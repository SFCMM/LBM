#ifndef LBM_BND_PERIODIC_H
#define LBM_BND_PERIODIC_H
#include "bnd_interface.h"
#include "constants.h"
#include "equilibrium_func.h"
#include "variables.h"

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBndCell_periodic : public LBMBndCell<LBTYPE> {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<LBEquation::Navier_Stokes, NDIM>;

 public:
  LBMBndCell_periodic() = delete;
  LBMBndCell_periodic(const GInt cellId, const VectorD<dim(LBTYPE)>& _normal) : LBMBndCell<LBTYPE>(cellId, _normal) {}
  virtual ~LBMBndCell_periodic() = default;

  // deleted constructors not needed
  //  LBMBndCell_periodic(const LBMBndCell_periodic&) = delete;
  //  LBMBndCell_periodic(LBMBndCell_periodic&&)      = delete;
  //  auto operator=(const LBMBndCell_periodic&) -> LBMBndCell_periodic& = delete;
  //  auto operator=(LBMBndCell_periodic&&) -> LBMBndCell_periodic& = delete;


  // todo: move to some common place
  void init(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surfConnected, GDouble pressure) {
    LBMBndCell<LBTYPE>::init();

    m_pressure               = pressure;
    const GDouble cellLength = surfConnected->cellLength(mapped());
    const auto    center     = surfConnected->center(mapped());
    const auto&   bb         = surfConnected->grid().boundingBox();

    // set which dists are to be set by this boundary condition
    for(GInt dist = 0; dist < NDIST; ++dist) {
      // determine if the dist points into the outside normal direction of bndry
      if(inDirection<dim(LBTYPE)>(normal(), LBMethod<LBTYPE>::m_dirs[dist])) {
        GBool inside = true;
        for(GInt dir = 0; dir < dim(LBTYPE); ++dir) {
          if(std::abs(normal()[dir]) > 0) {
            continue;
          }

          const GDouble coord = center[dir] + LBMethod<LBTYPE>::m_dirs[dist][dir] * cellLength;
          // todo: move to function isInside...
          if(coord < bb.min(dir) || coord > bb.max(dir)) {
            inside = false;
            break;
          }
        }

        if(inside) {
          m_bndIndex[m_noSetDists] = dist;
          ++m_noSetDists;
        }
      }
    }

    // todo: this needs to be moved since we cannot guarantee double matching cells are handled correctly
    for(GInt id = 0; id < m_noSetDists; ++id) {
      m_linkedCell[id]   = -1;
      const GInt dist    = m_bndIndex[id];
      auto       centerA = center;
      for(GInt dir = 0; dir < dim(LBTYPE); ++dir) {
        if(std::abs(normal()[dir]) > 0) {
          continue;
        }
        centerA[dir] += LBMethod<LBTYPE>::m_dirs[dist][dir] * cellLength;
      }

      //      GDouble minDist = std::numeric_limits<GDouble>::max();
      // connect cell of the boundary surface and connected surface
      for(const GInt cellIdB : surfConnected->getCellList()) {
        const auto& centerB = surfConnected->center(cellIdB);
        for(GInt dir = 0; dir < dim(LBTYPE); ++dir) {
          const GDouble diff = std::abs(centerA[dir] - centerB[dir]);
          if(diff <= GDoubleEps) {
            m_linkedCell[id] = cellIdB;
            break;
          }
        }
        if(m_linkedCell[id] >= 0) {
          break;
        }
      }

      ASSERT(m_linkedCell[id] >= 0, "Invalid Cellid");
    }

    for(GInt idA = 0; idA < m_noSetDists; ++idA) {
      for(GInt idB = idA + 1; idB < m_noSetDists; ++idB) {
        if(m_linkedCell[idA] == m_linkedCell[idB]) {
          TERMM(-1, "Invalid periodic bnd");
        }
      }
    }
  }

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) {
    if(!std::isnan(m_pressure)) {
      // need to calculate a pressure drop across the periodic boundary
      GDouble vsq = 0;
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vsq += vars(mapped(), VAR::velocity(dir)) * vars(mapped(), VAR::velocity(dir));
      }

      for(GInt id = 0; id < m_noSetDists; ++id) {
        const GInt dist = m_bndIndex[id];
        GDouble    cu   = vars(mapped(), VAR::velocity(0)) * LBMethod<LBTYPE>::m_dirs[dist][0];

        //        GDouble cu = 0;
        //        for(GInt dir = 0; dir < NDIM; ++dir) {
        //          cu += vars(mapped(), VAR::velocity(dir))  * LBMethod<LBTYPE>::m_dirs[dist][dir];
        //        }
        // for incompressible NSE rho = pressure
        fold(m_linkedCell[id], dist) =
            eq::defaultEq(LBMethod<LBTYPE>::m_weights[dist], m_pressure, cu, vsq) + f(mapped(), dist) - feq(mapped(), dist);
      }

    } else {
      //    do the same as in the propagation step
      for(GInt id = 0; id < m_noSetDists; ++id) {
        const GInt dist              = m_bndIndex[id];
        fold(m_linkedCell[id], dist) = f(mapped(), dist);
      }
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
             const std::function<GDouble&(GInt, GInt)>& /*vars*/) {}

 private:
  using LBMBndCell<LBTYPE>::mapped;
  using LBMBndCell<LBTYPE>::normal;
  using LBMBndCell<LBTYPE>::init;

  std::array<GInt, noDists(LBTYPE)> m_bndIndex{};
  std::array<GInt, noDists(LBTYPE)> m_linkedCell{};
  GInt                              m_noSetDists = 0;
  GDouble                           m_pressure   = NAN;
};

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBnd_Periodic : public LBMBndInterface {
 public:
  LBMBnd_Periodic(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surfConnected,
                  const json& properties) {
    for(const GInt cellId : surf->getCellList()) {
      m_bndCells.emplace_back(cellId, surf->normal(cellId));
    }
    if(config::has_config_value(properties, "pressure")) {
      logger << "Setting pressure for peridodic boundary" << std::endl;
      LBMBnd_Periodic<DEBUG_LEVEL, LBTYPE>::init(surfConnected, config::required_config_value<GDouble>(properties, "pressure"));
    } else {
      LBMBnd_Periodic<DEBUG_LEVEL, LBTYPE>::init(surfConnected);
    }
  }

  void init(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surfConnected, GDouble pressure = NAN) {
    for(auto& bndCell : m_bndCells) {
      bndCell.init(surfConnected, pressure);
    }
  }

  void initCnd(const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) override {
    // apply to all boundary cells
    for(auto& bndCell : m_bndCells) {
      bndCell.preApply(f, fold, feq, vars);
    }
  }
  void apply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
             const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

 private:
  std::vector<LBMBndCell_periodic<DEBUG_LEVEL, LBTYPE>> m_bndCells;
};
#endif // LBM_BND_PERIODIC_H
