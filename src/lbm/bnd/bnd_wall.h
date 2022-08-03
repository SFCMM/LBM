#ifndef LBM_BND_WALL_H
#define LBM_BND_WALL_H
#include "bnd_interface.h"
#include "bnd_wetnode.h"
#include "lbm/constants.h"
#include "lbm/moments.h"

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ, GBool TANGENTIALVELO>
class LBMBnd_wallBB : public LBMBnd_DirichletBB<DEBUG_LEVEL, LBTYPE, EQ> {
 private:
  static constexpr GInt NDIM = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NVAR = noVars<LBTYPE>(EQ);

 public:
  LBMBnd_wallBB(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties)
    : LBMBnd_DirichletBB<DEBUG_LEVEL, LBTYPE, EQ>(surf, properties, 0), m_bnd(surf) {
    if constexpr(TANGENTIALVELO) {
      m_tangentialVelo = config::required_config_value<GDouble>(properties, "tangentialVelocity");
      logger << "Setting tangentialVelocity " << m_tangentialVelo << std::endl;
    }
    LBMBnd_wallBB<DEBUG_LEVEL, LBTYPE, EQ, TANGENTIALVELO>::init();
  }
  ~LBMBnd_wallBB() override = default;

  // deleted constructors not needed
  LBMBnd_wallBB(const LBMBnd_wallBB&)                    = delete;
  LBMBnd_wallBB(LBMBnd_wallBB&&)                         = delete;
  auto operator=(const LBMBnd_wallBB&) -> LBMBnd_wallBB& = delete;
  auto operator=(LBMBnd_wallBB&&) -> LBMBnd_wallBB&      = delete;

  void init() {
    for(const GInt cellId : m_bnd->getCellList()) {
      m_dirichletValue.emplace_back();
      if constexpr(TANGENTIALVELO) {
        static constexpr GDouble parallelRequirement = 10 * GDoubleEps;


        const auto& _normal = VectorD<NDIM>(m_bnd->normal_p(cellId));
        for(GInt id = 0; id < noDists(LBTYPE); ++id) {
          // determine if the dist points into the outside normal direction of bndry
          if(inDirection<dim(LBTYPE)>(_normal, LBMethod<LBTYPE>::m_dirs[id])) {
            const GInt           insideDir = LBMethod<LBTYPE>::oppositeDist(id);
            VectorD<dim(LBTYPE)> dir;
            for(GInt axis = 0; axis < dim(LBTYPE); ++axis) {
              dir[axis] = LBMethod<LBTYPE>::m_dirs[insideDir][axis];
            }

            VectorD<dim(LBTYPE)> tangentialVec;
            if(dim(LBTYPE) == 2) {
              tangentialVec[0] = _normal[1];
              tangentialVec[1] = _normal[0];
            } else {
              TERMM(-1, "Not implemented");
            }


            const GDouble normalDotDir     = _normal.dot(dir);
            const GDouble tangentialDotDir = tangentialVec.dot(dir);
            const GDouble angle            = gcem::acos(normalDotDir / dir.norm());
            const GBool   parallel         = std::abs(angle - PI) < parallelRequirement;

            // vector is parallel to the normal vector -> velocity = 0
            if(!parallel) {
              m_dirichletValue.back()[insideDir] = m_tangentialVelo * tangentialDotDir;
            }
          }
        }
      } else {
        m_dirichletValue.back().fill(0);
      }
    }
  }

  void initCnd(const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) override {
    GInt index = 0;
    for(const GInt cellId : m_bnd->getCellList()) {
      // TANGENTIALVELO is false when the value is zero
      LBMBnd_DirichletBB<DEBUG_LEVEL, LBTYPE, EQ>::template apply<!TANGENTIALVELO, true>(cellId, &m_dirichletValue[index][0], fpre, fold,
                                                                                         feq, vars);
      ++index;
    }
  }

 private:
  GDouble                                                                               m_tangentialVelo = 0.1;
  std::vector<std::array<GDouble, noDists(LBTYPE) * static_cast<GInt>(TANGENTIALVELO)>> m_dirichletValue;

  const SurfaceInterface* m_bnd = nullptr;
};

/// Equilibrium wall boundary condition.
/// \tparam DEBUG_LEVEL
/// \tparam LBTYPE
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBnd_wallEq : public LBMBndInterface, protected LBMBnd_wallWetnode<DEBUG_LEVEL, LBTYPE> {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<LBEquationType::Navier_Stokes, NDIM>;

 public:
  LBMBnd_wallEq(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties)
    : LBMBnd_wallWetnode<DEBUG_LEVEL, LBTYPE>(surf), m_bnd(surf) {
    if(!config::has_config_value(properties, "velocity")) {
      m_apply = &LBMBnd_wallEq::apply_0;
      m_wallV.fill(0);
      logger << " Wall with no-slip condition" << std::endl;
    } else {
      m_apply = &LBMBnd_wallEq::apply_constV;
      m_wallV = config::required_config_value<NDIM>(properties, "velocity");
      logger << " Wall with constant velocity" << std::endl;
    }
  }
  LBMBnd_wallEq(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf) : LBMBnd_wallWetnode<DEBUG_LEVEL, LBTYPE>(surf), m_bnd(surf) {}
  ~LBMBnd_wallEq() override = default;

  // deleted constructors not needed
  LBMBnd_wallEq(const LBMBnd_wallEq&)                    = delete;
  LBMBnd_wallEq(LBMBnd_wallEq&&)                         = delete;
  auto operator=(const LBMBnd_wallEq&) -> LBMBnd_wallEq& = delete;
  auto operator=(LBMBnd_wallEq&&) -> LBMBnd_wallEq&      = delete;

  void initCnd(const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>&    fold,
             const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& vars) override {
    m_apply(this, fpre, fold, vars);
  }

 protected:
  virtual void apply_0(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& fold,
                       const std::function<GDouble&(GInt, GInt)>& vars) {
    set_wallV(vars);

    // update the density value with velocity set to 0
    GInt index = 0;
    for(const auto cellId : m_bnd->getCellList()) {
      calcDensity_limited<NDIM, NDIST, LBEquationType::Navier_Stokes, true>(cellId, this->limitedDist()[index], this->limitedConst()[index],
                                                                            nullptr, fold, vars);
      ++index;
    }

    // for wall with 0 velocity
    for(const auto cellId : m_bnd->getCellList()) {
      for(GInt dist = 0; dist < NDIST; ++dist) {
        // todo: make settable
        fold(cellId, dist) = eq::defaultEq(method::m_weights[dist], vars(cellId, VAR::rho()));
      }
    }
  }

  virtual void apply_constV(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& fold,
                            const std::function<GDouble&(GInt, GInt)>& vars) {
    set_wallV(vars);

    // update the density value with velocity set to the wall velocity
    GInt index = 0;
    for(const auto cellId : m_bnd->getCellList()) {
      calcDensity_limited<NDIM, NDIST, LBEquationType::Navier_Stokes, false>(
          cellId, this->limitedDist()[index], this->limitedConst()[index], &this->normal()[index][0], fold, vars);
      ++index;
    }

    // for wall with constant velocity
    for(const auto cellId : m_bnd->getCellList()) {
      // todo: make settable
      eq::defaultEq<LBTYPE>(&fold(cellId, 0), vars(cellId, VAR::rho()), &vars(cellId, VAR::velocity(0)));
    }
  }

  [[nodiscard]] auto bnd() const -> const SurfaceInterface* { return m_bnd; }
  auto               wallV() -> VectorD<NDIM>& { return m_wallV; }

 private:
  void set_wallV(const std::function<GDouble&(GInt, GInt)>& vars) {
    // at wall set wall velocity
    for(const auto cellId : m_bnd->getCellList()) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(cellId, VAR::velocity(dir)) = m_wallV[dir];
      }
    }
  }

  const SurfaceInterface* m_bnd = nullptr;

  std::function<void(LBMBnd_wallEq*, const std::function<GDouble&(GInt, GInt)>&, const std::function<GDouble&(GInt, GInt)>&,
                     const std::function<GDouble&(GInt, GInt)>&)>
      m_apply;

  VectorD<NDIM> m_wallV;
};

/// Non-equilibrium extrapolation wall boundary condition.
/// \tparam DEBUG_LEVEL
/// \tparam LBTYPE
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ>
class LBMBnd_wallNEEM : /*public LBMBndInterface,*/ public LBMBnd_wallEq<DEBUG_LEVEL, LBTYPE> {
 private:
  using METH                  = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<EQ, NDIM>;

 public:
  LBMBnd_wallNEEM(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties) : LBMBnd_wallEq<DEBUG_LEVEL, LBTYPE>(surf) {
    if(!config::has_config_value(properties, "velocity")) {
      m_apply = &LBMBnd_wallNEEM::apply_0NEEM;
      this->wallV().fill(0);
      logger << " NEEM with no-slip condition" << std::endl;
    } else {
      m_apply       = &LBMBnd_wallNEEM::apply_constVNEEM;
      this->wallV() = config::required_config_value<NDIM>(properties, "velocity");
      logger << " NEEM with constant velocity" << std::endl;
    }
    init();
  }
  ~LBMBnd_wallNEEM() override = default;

  // deleted constructors not needed
  LBMBnd_wallNEEM(const LBMBnd_wallNEEM&)                    = delete;
  LBMBnd_wallNEEM(LBMBnd_wallNEEM&&)                         = delete;
  auto operator=(const LBMBnd_wallNEEM&) -> LBMBnd_wallNEEM& = delete;
  auto operator=(LBMBnd_wallNEEM&&) -> LBMBnd_wallNEEM&      = delete;

  void init() {
    auto extrapolationDir = [&](const GDouble* normal) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        if(normal[dir] < 0) {
          return 2 * dir + 1;
        }
        if(normal[dir] > 0) {
          return 2 * dir;
        }
      }
      return static_cast<GInt>(-1);
    };

    for(const auto bndCellId : this->bnd()->getCellList()) {
      const GDouble* normal    = this->bnd()->normal_p(bndCellId);
      const GInt     extDir    = extrapolationDir(normal);
      const GInt     extCellId = this->bnd()->neighbor(bndCellId, extDir);

      // handle corner
      //      if(NDIM == 2) {
      //        if((this->bnd()->neighbor(bndCellId, 0) == INVALID_CELLID || this->bnd()->neighbor(bndCellId, 1) == INVALID_CELLID)
      //           && (this->bnd()->neighbor(bndCellId, 2) == INVALID_CELLID || this->bnd()->neighbor(bndCellId, 3) == INVALID_CELLID)) {
      //          // place an invalid cellId to mark corner
      //          m_extrapolationCellId.emplace_back(INVALID_CELLID);
      //          continue;
      //        }
      //      }
      cerr0 << "bndCellId " << bndCellId << " extrapolationDir " << extDir << " extracellId " << extCellId << std::endl;


      m_extrapolationCellId.emplace_back(extCellId);
      if(extCellId == INVALID_CELLID) {
        cerr0 << "bndCellId " << bndCellId << " extrapolationDir " << extDir << std::endl;
        TERMM(-1, "No valid extrapolation cellId");
      }
    }
  }

  void initCnd(const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) override {
    m_apply(this, fpre, fold, feq, vars);
  }

 private:
  void apply_0NEEM(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold,
                   const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) {
    LBMBnd_wallEq<DEBUG_LEVEL, LBTYPE>::apply_0(fpre, fold, vars);

    apply_general(fpre, fold, feq, vars);
  }

  void apply_constVNEEM(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold,
                        const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) {
    LBMBnd_wallEq<DEBUG_LEVEL, LBTYPE>::apply_constV(fpre, fold, vars);

    apply_general(fpre, fold, feq, vars);
  }

  void apply_general(const std::function<GDouble&(GInt, GInt)>& /*fpre*/, const std::function<GDouble&(GInt, GInt)>& fold,
                     const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>&  vars) {
    if constexpr(EQ == LBEquationType::Navier_Stokes) {
      // recalculate density and velocity dist because fold is post streaming!
      calcDensity<NDIM, NDIST, EQ>(m_extrapolationCellId, fold, vars);
      calcVelocity<NDIM, NDIST, EQ>(m_extrapolationCellId, fold, vars);

      // add the non-equilibrium part from the extrapolation cell
      GInt index = 0;
      for(const auto cellId : this->bnd()->getCellList()) {
        const GInt extCellId = m_extrapolationCellId[index];

        //  just set equilibrium dist if there is no valid extrapolation cell
        if(extCellId != INVALID_CELLID) {
          for(GInt dist = 0; dist < NDIST; ++dist) {
            fold(cellId, dist) +=
                fold(extCellId, dist) - eq::defaultEq<LBTYPE>(dist, vars(extCellId, VAR::rho()), &vars(extCellId, VAR::velocity(0)));
          }
        }
        ++index;
      }
    }
  }

  std::function<void(LBMBnd_wallNEEM*, const std::function<GDouble&(GInt, GInt)>&, const std::function<GDouble&(GInt, GInt)>&,
                     const std::function<GDouble&(GInt, GInt)>&, const std::function<GDouble&(GInt, GInt)>&)>
      m_apply;

  std::vector<GInt> m_extrapolationCellId;
};


// todo: add literature info
/// Nonequilibrium bounce-back model can be *only* used for D2Q9. In the literature there is doubt that it is useful in 3D.
/// The model leads to a relaxation independent result!
/// \tparam DEBUG_LEVEL
/// \tparam LBTYPE
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBnd_wallNEBB : public LBMBndInterface, public LBMBnd_wallWetnode<DEBUG_LEVEL, LBTYPE> {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<LBEquationType::Navier_Stokes, NDIM>;

 public:
  LBMBnd_wallNEBB(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties)
    : LBMBnd_wallWetnode<DEBUG_LEVEL, LBTYPE>(surf), m_bnd(surf) {
    if(LBTYPE != LBMethodType::D2Q9) {
      // Note: there is doubt in the current literature that this type of boundary condition make sense in 3D
      TERMM(-1, "Not implemented for this distribution!");
    }

    if(!config::has_config_value(properties, "velocity")) {
      m_apply = &LBMBnd_wallNEBB::apply_0;
      logger << " NEBB with no-slip condition" << std::endl;
    } else {
      m_apply = &LBMBnd_wallNEBB::apply_constV;
      m_wallV = config::required_config_value<NDIM>(properties, "velocity");
      logger << " NEBB with constant velocity" << std::endl;
    }
  }
  ~LBMBnd_wallNEBB() override = default;

  // deleted constructors not needed
  LBMBnd_wallNEBB(const LBMBnd_wallNEBB&)                    = delete;
  LBMBnd_wallNEBB(LBMBnd_wallNEBB&&)                         = delete;
  auto operator=(const LBMBnd_wallNEBB&) -> LBMBnd_wallNEBB& = delete;
  auto operator=(LBMBnd_wallNEBB&&) -> LBMBnd_wallNEBB&      = delete;

  void initCnd(const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) override {
    m_apply(this, fpre, fold, feq, vars);
  }


 private:
  void apply_0(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>&   fold,
               const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& vars) {
    // at wall set 0 velocity
    for(const auto cellId : m_bnd->getCellList()) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(cellId, VAR::velocity(dir)) = 0;
      }
    }

    GInt index = 0;
    for(const auto cellId : m_bnd->getCellList()) {
      calcDensity_limited<NDIM, NDIST, LBEquationType::Navier_Stokes, true>(cellId, this->limitedDist()[index], this->limitedConst()[index],
                                                                            nullptr, fold, vars);
      ++index;
    }

    // reflect main directions
    for(const auto bndCellId : m_bnd->getCellList()) {
      for(GInt dist = 0; dist < cartesian::maxNoNghbrs<NDIM>(); ++dist) {
        if(m_bnd->neighbor(bndCellId, dist) == INVALID_CELLID) {
          fold(bndCellId, cartesian::oppositeDir(dist)) = fold(bndCellId, dist);
        }
      }
    }

    // reflect diagonal direction
    index = 0;
    for(const auto bndCellId : m_bnd->getCellList()) {
      const auto& normal = this->normal()[index];
      if(normal[0] < 0) {
        // +x
        fold(bndCellId, 6) = fold(bndCellId, 4) + 0.5 * (fold(bndCellId, 3) - fold(bndCellId, 2));
        fold(bndCellId, 7) = fold(bndCellId, 5) - 0.5 * (fold(bndCellId, 3) - fold(bndCellId, 2));
      } else if(normal[0] > 0) {
        // -x
        fold(bndCellId, 4) = fold(bndCellId, 6) - 0.5 * (fold(bndCellId, 3) - fold(bndCellId, 2));
        fold(bndCellId, 5) = fold(bndCellId, 7) + 0.5 * (fold(bndCellId, 3) - fold(bndCellId, 2));
      } else if(normal[1] < 0) {
        // +y
        fold(bndCellId, 5) = fold(bndCellId, 7) - 0.5 * (fold(bndCellId, 1) - fold(bndCellId, 0));
        fold(bndCellId, 6) = fold(bndCellId, 4) + 0.5 * (fold(bndCellId, 1) - fold(bndCellId, 0));
      } else if(normal[1] > 0) {
        //-y
        fold(bndCellId, 7) = fold(bndCellId, 5) + 0.5 * (fold(bndCellId, 1) - fold(bndCellId, 0));
        fold(bndCellId, 4) = fold(bndCellId, 6) - 0.5 * (fold(bndCellId, 1) - fold(bndCellId, 0));
      }
    }
    ++index;
  }

  void apply_constV(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>&   fold,
                    const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& vars) {
    // todo: this is not correct see olb-1.4/src/boundary/wallFunctionB*/computeFneqRNEBB
    GInt index = 0;
    for(const auto cellId : m_bnd->getCellList()) {
      calcDensity_limited<NDIM, NDIST, LBEquationType::Navier_Stokes, false>(
          cellId, this->limitedDist()[index], this->limitedConst()[index], &this->normal()[index][0], fold, vars);
      ++index;
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(cellId, VAR::velocity(dir)) = m_wallV[dir];
      }
    }

    index = 0;
    for(const auto bndCellId : m_bnd->getCellList()) {
      const auto& normal = this->normal()[index];
      // reflect main direction
      for(GInt dist = 0; dist < cartesian::maxNoNghbrs<NDIM>(); ++dist) {
        const GInt dir = dist / 2;
        if(m_bnd->neighbor(bndCellId, dist) == INVALID_CELLID && std::abs(normal[dir]) > 0) {
          fold(bndCellId, cartesian::oppositeDir(dist)) = fold(bndCellId, dist);
          if(normal[dir] < 0) {
            fold(bndCellId, cartesian::oppositeDir(dist)) -= 2.0 / 3.0 * vars(bndCellId, VAR::rho()) * m_wallV[dir];
          } else {
            fold(bndCellId, cartesian::oppositeDir(dist)) += 2.0 / 3.0 * vars(bndCellId, VAR::rho()) * m_wallV[dir];
          }
        }
      }
      ++index;
    }

    // reflect diagonal direction
    index = 0;
    for(const auto bndCellId : m_bnd->getCellList()) {
      const GDouble rho    = vars(bndCellId, VAR::rho());
      const auto&   normal = this->normal()[index];
      if(normal[0] > 0) {
        // +x right wall
        fold(bndCellId, 6) =
            fold(bndCellId, 4) + 0.5 * (fold(bndCellId, 3) - fold(bndCellId, 2)) - 0.5 * rho * m_wallV[1] - 1.0 / 6.0 * rho * m_wallV[0];
        fold(bndCellId, 7) =
            fold(bndCellId, 5) - 0.5 * (fold(bndCellId, 3) - fold(bndCellId, 2)) + 0.5 * rho * m_wallV[1] - 1.0 / 6.0 * rho * m_wallV[0];
      } else if(normal[0] < 0) {
        // -x left wall
        fold(bndCellId, 4) =
            fold(bndCellId, 6) - 0.5 * (fold(bndCellId, 3) - fold(bndCellId, 2)) + 0.5 * rho * m_wallV[1] + 1.0 / 6.0 * rho * m_wallV[0];
        fold(bndCellId, 5) =
            fold(bndCellId, 7) + 0.5 * (fold(bndCellId, 3) - fold(bndCellId, 2)) - 0.5 * rho * m_wallV[1] + 1.0 / 6.0 * rho * m_wallV[0];
      } else if(normal[1] > 0) {
        // +y top wall
        fold(bndCellId, 5) =
            fold(bndCellId, 7) - 0.5 * (fold(bndCellId, 1) - fold(bndCellId, 0)) + 0.5 * rho * m_wallV[0] - 1.0 / 6.0 * rho * m_wallV[1];
        fold(bndCellId, 6) =
            fold(bndCellId, 4) + 0.5 * (fold(bndCellId, 1) - fold(bndCellId, 0)) - 0.5 * rho * m_wallV[0] - 1.0 / 6.0 * rho * m_wallV[1];
      } else if(normal[1] < 0) {
        //-y bottom wall
        fold(bndCellId, 7) =
            fold(bndCellId, 5) + 0.5 * (fold(bndCellId, 1) - fold(bndCellId, 0)) - 0.5 * rho * m_wallV[0] + 1.0 / 6.0 * rho * m_wallV[1];
        fold(bndCellId, 4) =
            fold(bndCellId, 6) - 0.5 * (fold(bndCellId, 1) - fold(bndCellId, 0)) + 0.5 * rho * m_wallV[0] + 1.0 / 6.0 * rho * m_wallV[1];
      }
      ++index;
    }
  }

  std::function<void(LBMBnd_wallNEBB*, const std::function<GDouble&(GInt, GInt)>&, const std::function<GDouble&(GInt, GInt)>&,
                     const std::function<GDouble&(GInt, GInt)>&, const std::function<GDouble&(GInt, GInt)>&)>
      m_apply;

  const SurfaceInterface* m_bnd = nullptr;
  VectorD<NDIM>           m_wallV;
};
#endif // LBM_BND_WALL_H
