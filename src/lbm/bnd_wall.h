#ifndef LBM_BND_WALL_H
#define LBM_BND_WALL_H
#include "bnd_interface.h"
#include "constants.h"
#include "moments.h"

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

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) {
    // iterate over the distributions that need to be set
    for(GInt id = 0; id < m_noSetDists; ++id) {
      // dists that point into the outside direction
      const GInt dist = m_bndIndex[id];
      // dist in reflected/opposite direction i.e. to the inside of the wall
      GInt oppositeDist = LBMethod<LBTYPE>::oppositeDist(dist);

      // standard bounceback i.e. distribution that hits the wall is reflected to the opposite distribution direction
      fold(mapped(), oppositeDist) = f(mapped(), dist);

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
        for(GInt n = 0; n < dim(LBTYPE); ++n) {
          dir[n] = LBMethod<LBTYPE>::m_dirs[dist][n];
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

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {
    for(auto& bndCell : m_bndCells) {
      bndCell.apply(f, fold);
    }
  }

 private:
  GDouble                                                m_tangentialVelo = 0.1;
  std::vector<LBMBndCell_wallBB<LBTYPE, TANGENTIALVELO>> m_bndCells;
};

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBnd_wallEq : public LBMBndInterface {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<LBEquation::Navier_Stokes, NDIM>;

 public:
  LBMBnd_wallEq(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties) : m_bnd(surf) {
    if(!config::has_config_value(properties, "velocity")) {
      m_apply = &LBMBnd_wallEq::apply_0;
      logger << " Wall with no-slip condition" << std::endl;
    } else {
      m_apply       = &LBMBnd_wallEq::apply_constV;
      m_tangentialV = config::required_config_value<NDIM>(properties, "velocity");
      logger << " Wall with constant velocity" << std::endl;
    }
  }
  ~LBMBnd_wallEq() override = default;

  // deleted constructors not needed
  LBMBnd_wallEq(const LBMBnd_wallEq&) = delete;
  LBMBnd_wallEq(LBMBnd_wallEq&&)      = delete;
  auto operator=(const LBMBnd_wallEq&) -> LBMBnd_wallEq& = delete;
  auto operator=(LBMBnd_wallEq&&) -> LBMBnd_wallEq& = delete;

  void initCnd(const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& vars) override {
    m_apply(this, f, fold, vars);
  }

 private:
  void apply_0(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& fold,
               const std::function<GDouble&(GInt, GInt)>& vars) {
    //     at wall set 0 velocity
    for(const auto cellId : m_bnd->getCellList()) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(cellId, VAR::velocity(dir)) = 0;
      }
    }

    for(const auto cellId : m_bnd->getCellList()) {
      vars(cellId, VAR::rho()) =
          1.0 / (1.0 - vars(cellId, VAR::velocity(1)))
          * (fold(cellId, 8) + fold(cellId, 1) + fold(cellId, 0) + 2.0 * (fold(cellId, 6) + fold(cellId, 2) + fold(cellId, 5)));
    }

    // for wall with 0 velocity
    for(const auto cellId : m_bnd->getCellList()) {
      for(GInt dist = 0; dist < NDIST; ++dist) {
        // todo: make settable
        fold(cellId, dist) = eq::defaultEq(method::m_weights[dist], vars(cellId, VAR::rho()), 0, 0);
      }
    }
  }

  void apply_constV(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& fold,
                    const std::function<GDouble&(GInt, GInt)>& vars) {
    // update the density in the boundary cells
    //    calcDensity<NDIM, NDIST, LBEquation::Navier_Stokes>(m_bnd->getCellList(), fold, vars);


    for(const auto cellId : m_bnd->getCellList()) {
      vars(cellId, VAR::rho()) =
          1.0 / (1.0 - vars(cellId, VAR::velocity(1)))
          * (fold(cellId, 8) + fold(cellId, 1) + fold(cellId, 0) + 2.0 * (fold(cellId, 6) + fold(cellId, 2) + fold(cellId, 5)));
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(cellId, VAR::velocity(dir)) = m_tangentialV[dir];
      }
    }

    // for wall with 0 velocity
    for(const auto cellId : m_bnd->getCellList()) {
      GDouble vsq = 0;
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vsq += vars(cellId, VAR::velocity(dir)) * vars(cellId, VAR::velocity(dir));
      }

      for(GInt dist = 0; dist < NDIST; ++dist) {
        GDouble cu = 0;
        for(GInt dir = 0; dir < NDIM; ++dir) {
          cu += vars(cellId, VAR::velocity(dir)) * LBMethod<LBTYPE>::m_dirs[dist][dir];
        }
        // todo: make settable
        fold(cellId, dist) = eq::defaultEq(method::m_weights[dist], vars(cellId, VAR::rho()), cu, vsq);
      }
    }
  }

  const SurfaceInterface* m_bnd = nullptr;
  std::function<void(LBMBnd_wallEq*, const std::function<GDouble&(GInt, GInt)>&, const std::function<GDouble&(GInt, GInt)>&,
                     const std::function<GDouble&(GInt, GInt)>&)>
      m_apply;

  VectorD<NDIM> m_tangentialV;
};
#endif // LBM_BND_WALL_H
