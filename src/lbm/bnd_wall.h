#ifndef LBM_BND_WALL_H
#define LBM_BND_WALL_H
#include "bnd_interface.h"
#include "constants.h"
#include "moments.h"

// todo: remove!
template <LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBndCell_wallBB : public LBMBndCell<LBTYPE> {
  // class LBMBndCell_wallBB {
 public:
  LBMBndCell_wallBB() = delete;
  //  LBMBndCell_wallBB() = default;
  LBMBndCell_wallBB(const GInt cellId, const VectorD<dim(LBTYPE)>& _normal) : LBMBndCell<LBTYPE>(cellId, _normal) {}
  //  LBMBndCell_wallBB(const GInt cellId) {}
  ~LBMBndCell_wallBB() override = default;

  void init() override {
    LBMBndCell<LBTYPE>::init();

    // precalculate weight and bndIndex
    for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
      // determine if the dist points into the outside normal direction of bndry
      if(inDirection<dim(LBTYPE)>(normal(), LBMethod<LBTYPE>::m_dirs[dist])) {
        m_bndIndex[m_noSetDists] = dist;
        ++m_noSetDists;
      }
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold) {
    // iterate over the distributions that need to be set
    for(GInt id = 0; id < m_noSetDists; ++id) {
      // dists that point into the outside direction
      const GInt dist = m_bndIndex[id];
      // dist in reflected/opposite direction i.e. to the inside of the wall
      GInt oppositeDist = LBMethod<LBTYPE>::oppositeDist(dist);

      // standard bounceback i.e. distribution that hits the wall is reflected to the opposite distribution direction
      fold(mapped(), oppositeDist) = fpre(mapped(), dist);

      if constexpr(TANGENTIALVELO) {
        // todo: actually calculate this
        GDouble avg_density = 1.0;
        fold(mapped(), oppositeDist) += avg_density * m_tangentialVelo[oppositeDist];
      }
    }
  }

  void setTangentialVelocity(const GDouble tangentialVelo) {
    if constexpr(TANGENTIALVELO) {
      static constexpr GDouble parallelRequirement = 10 * GDoubleEps;

      m_tangentialVelo.fill(0);
      const auto& _normal = normal();
      for(GInt id = 0; id < m_noSetDists; ++id) {
        const GInt           dist = LBMethod<LBTYPE>::oppositeDist(m_bndIndex[id]);
        VectorD<dim(LBTYPE)> dir;
        for(GInt axis = 0; axis < dim(LBTYPE); ++axis) {
          dir[axis] = LBMethod<LBTYPE>::m_dirs[dist][axis];
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
          m_tangentialVelo[dist] = tangentialDotDir * 2 * LBMethod<LBTYPE>::m_weights[dist] * 1.0 / lbm_cssq * tangentialVelo;
        }
      }
    }
  }

  // todo: fix me
  //  deleted constructors not needed
  //  LBMBndCell_wallBB(const LBMBndCell_wallBB&) = delete;
  //  LBMBndCell_wallBB(LBMBndCell_wallBB&&)      = delete;
  //  auto operator=(const LBMBndCell_wallBB&) -> LBMBndCell_wallBB& = delete;
  //  auto operator=(LBMBndCell_wallBB&&) -> LBMBndCell_wallBB& = delete;

 private:
  using LBMBndCell<LBTYPE>::mapped;
  using LBMBndCell<LBTYPE>::normal;

  std::array<GDouble, noDists(LBTYPE) * static_cast<GInt>(TANGENTIALVELO)> m_tangentialVelo;
  std::array<GInt, noDists(LBTYPE)>                                        m_bndIndex;
  GInt                                                                     m_noSetDists = 0;
};

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBnd_wallBB : public LBMBndInterface {
 public:
  LBMBnd_wallBB(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties) {
    for(const GInt cellId : surf->getCellList()) {
      m_bndCells.emplace_back(cellId, surf->normal(cellId));
    }
    if constexpr(TANGENTIALVELO) {
      m_tangentialVelo = config::required_config_value<GDouble>(properties, "tangentialVelocity");
      logger << "Setting tangentialVelocity " << m_tangentialVelo << std::endl;
    }
    LBMBnd_wallBB<DEBUG_LEVEL, LBTYPE, TANGENTIALVELO>::init();
  }
  ~LBMBnd_wallBB() override = default;

  // deleted constructors not needed
  LBMBnd_wallBB(const LBMBnd_wallBB&) = delete;
  LBMBnd_wallBB(LBMBnd_wallBB&&)      = delete;
  auto operator=(const LBMBnd_wallBB&) -> LBMBnd_wallBB& = delete;
  auto operator=(LBMBnd_wallBB&&) -> LBMBnd_wallBB& = delete;

  void init() {
    for(auto& bndCell : m_bndCells) {
      bndCell.init();
      bndCell.setTangentialVelocity(m_tangentialVelo);
    }
  }

  void initCnd(const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {
    for(auto& bndCell : m_bndCells) {
      bndCell.apply(fpre, fold);
    }
  }

 private:
  GDouble                                                m_tangentialVelo = 0.1;
  std::vector<LBMBndCell_wallBB<LBTYPE, TANGENTIALVELO>> m_bndCells;
};

/// Base class for wet node boundary conditions.
/// \tparam DEBUG_LEVEL
/// \tparam LBTYPE
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBnd_wallWetnode {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<LBEquationType::Navier_Stokes, NDIM>;

 public:
  LBMBnd_wallWetnode(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf) {
    // determine limited dist set etc.
    for(const GInt cellId : surf->getCellList()) {
      m_limDist.emplace_back();
      m_limConst.emplace_back();
      m_limConst.back().fill(0);
      m_normal.emplace_back(VectorD<NDIM>(surf->normal_p(cellId)));
      for(GInt dir = 0; dir < NDIST - 1; ++dir) {
        const GInt oppositeDist = LBMethod<LBTYPE>::oppositeDist(dir);

        if(inDirection<dim(LBTYPE)>(m_normal.back(), LBMethod<LBTYPE>::m_dirs[dir])
           && surf->neighbor(cellId, oppositeDist) != INVALID_CELLID) {
          m_limDist.back().emplace(dir);

          m_limConst.back()[dir] = 2;

        } else if(orthogonal<dim(LBTYPE)>(m_normal.back(), LBMethod<LBTYPE>::m_dirs[dir])
                  && surf->neighbor(cellId, oppositeDist) != INVALID_CELLID) {
          m_limDist.back().emplace(dir);

          // cell is a corner
          if(surf->neighbor(cellId, dir) == INVALID_CELLID) {
            // this is only correct for wall velocity = 0
            m_limConst.back()[dir] = 2;
          } else {
            m_limConst.back()[dir] = 1;
          }
        }
      }
      const GInt sumC = std::accumulate(m_limConst.back().begin(), m_limConst.back().end(), 0);
      if(sumC != NDIST - 1) {
        // we clear to mark a corner
        m_limDist.back().clear();
      } else {
        m_limDist.back().emplace(NDIST - 1);
        m_limConst.back()[NDIST - 1] = 1;
      }
    }
  }

  [[nodiscard]] auto limitedDist() const -> const std::vector<std::set<GInt>>& { return m_limDist; }

  auto limitedConst() const -> const std::vector<std::array<GDouble, NDIST>>& { return m_limConst; }

  auto normal() const -> const std::vector<VectorD<NDIM>>& { return m_normal; }

 private:
  std::vector<std::set<GInt>>             m_limDist;
  std::vector<std::array<GDouble, NDIST>> m_limConst;
  std::vector<VectorD<NDIM>>              m_normal;
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
  LBMBnd_wallEq(const LBMBnd_wallEq&) = delete;
  LBMBnd_wallEq(LBMBnd_wallEq&&)      = delete;
  auto operator=(const LBMBnd_wallEq&) -> LBMBnd_wallEq& = delete;
  auto operator=(LBMBnd_wallEq&&) -> LBMBnd_wallEq& = delete;

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
    // at wall set 0 velocity
    for(const auto cellId : m_bnd->getCellList()) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(cellId, VAR::velocity(dir)) = 0;
      }
    }

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
        fold(cellId, dist) = eq::defaultEq(method::m_weights[dist], vars(cellId, VAR::rho()), 0, 0);
      }
    }
  }

  virtual void apply_constV(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& fold,
                            const std::function<GDouble&(GInt, GInt)>& vars) {
    GInt index = 0;
    for(const auto cellId : m_bnd->getCellList()) {
      calcDensity_limited<NDIM, NDIST, LBEquationType::Navier_Stokes, false>(
          cellId, this->limitedDist()[index], this->limitedConst()[index], &this->normal()[index][0], fold, vars);
      ++index;
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(cellId, VAR::velocity(dir)) = m_wallV[dir];
      }
    }

    // for wall with constant velocity
    for(const auto cellId : m_bnd->getCellList()) {
      // todo: make settable
      eq::defaultEq<LBTYPE>(&fold(cellId, 0), vars(cellId, VAR::rho()), &vars(cellId, VAR::velocity(0)));
    }
  }

  [[nodiscard]] auto bnd() const -> const SurfaceInterface* { return m_bnd; }
  auto               v_wall() -> VectorD<NDIM>& { return m_wallV; }

 private:
  const SurfaceInterface* m_bnd = nullptr;

  std::function<void(LBMBnd_wallEq*, const std::function<GDouble&(GInt, GInt)>&, const std::function<GDouble&(GInt, GInt)>&,
                     const std::function<GDouble&(GInt, GInt)>&)>
      m_apply;

  VectorD<NDIM> m_wallV;
};

/// Non-equilibrium extrapolation wall boundary condition.
/// \tparam DEBUG_LEVEL
/// \tparam LBTYPE
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBnd_wallNEEM : /*public LBMBndInterface,*/ public LBMBnd_wallEq<DEBUG_LEVEL, LBTYPE> {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<LBEquationType::Navier_Stokes, NDIM>;

 public:
  LBMBnd_wallNEEM(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties) : LBMBnd_wallEq<DEBUG_LEVEL, LBTYPE>(surf) {
    if(!config::has_config_value(properties, "velocity")) {
      m_apply = &LBMBnd_wallNEEM::apply_0NEEM;
      logger << " NEEM with no-slip condition" << std::endl;
    } else {
      m_apply        = &LBMBnd_wallNEEM::apply_constVNEEM;
      this->v_wall() = config::required_config_value<NDIM>(properties, "velocity");
      logger << " NEEM with constant velocity" << std::endl;
    }
    init();
  }
  ~LBMBnd_wallNEEM() override = default;

  // deleted constructors not needed
  LBMBnd_wallNEEM(const LBMBnd_wallNEEM&) = delete;
  LBMBnd_wallNEEM(LBMBnd_wallNEEM&&)      = delete;
  auto operator=(const LBMBnd_wallNEEM&) -> LBMBnd_wallNEEM& = delete;
  auto operator=(LBMBnd_wallNEEM&&) -> LBMBnd_wallNEEM& = delete;

  void init() {
    for(const auto bndCellId : this->bnd()->getCellList()) {
      GInt extrapolationDir = -1;
      for(GInt dist = 0; dist < cartesian::maxNoNghbrs<NDIM>(); ++dist) {
        if(this->bnd()->neighbor(bndCellId, dist) == INVALID_CELLID) {
          if(extrapolationDir < 0) {
            extrapolationDir = dist;
          } else {
            // todo: move to function
            //  we are at a corner so we have to use a diagonal neighbor
            if(extrapolationDir == 0 && dist == 2) {
              extrapolationDir = 6;
            }
            if(extrapolationDir == 0 && dist == 3) {
              extrapolationDir = 7;
            }
            if(extrapolationDir == 1 && dist == 3) {
              extrapolationDir = 4;
            }
            if(extrapolationDir == 1 && dist == 2) {
              extrapolationDir = 5;
            }
          }
        }
      }
      const GInt extrapolationCellId = this->bnd()->neighbor(bndCellId, cartesian::oppositeDir<NDIM>(extrapolationDir));
      m_extrapolationCellId.emplace_back(extrapolationCellId);
      if(extrapolationCellId == INVALID_CELLID) {
        cerr0 << "bndCellId " << bndCellId << " extrapolationDir " << extrapolationDir << std::endl;
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

    // for wall with 0 velocity
    GInt index = 0;
    for(const auto cellId : this->bnd()->getCellList()) {
      const GInt extraPolationCellId = m_extrapolationCellId[index];
      for(GInt dist = 0; dist < NDIST; ++dist) {
        // todo: make settable
        fold(cellId, dist) += fold(extraPolationCellId, dist) - feq(extraPolationCellId, dist);
      }
      ++index;
    }
  }

  void apply_constVNEEM(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold,
                        const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) {
    LBMBnd_wallEq<DEBUG_LEVEL, LBTYPE>::apply_constV(fpre, fold, vars);

    // for wall with constant velocity
    GInt index = 0;
    for(const auto cellId : this->bnd()->getCellList()) {
      const GInt extraPolationCellId = m_extrapolationCellId[index];
      for(GInt dist = 0; dist < NDIST; ++dist) {
        fold(cellId, dist) += fold(extraPolationCellId, dist) - feq(extraPolationCellId, dist);
      }
      ++index;
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
  LBMBnd_wallNEBB(const LBMBnd_wallNEBB&) = delete;
  LBMBnd_wallNEBB(LBMBnd_wallNEBB&&)      = delete;
  auto operator=(const LBMBnd_wallNEBB&) -> LBMBnd_wallNEBB& = delete;
  auto operator=(LBMBnd_wallNEBB&&) -> LBMBnd_wallNEBB& = delete;

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

    // reflect main direction
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
