#ifndef LBM_BND_PERIODIC_H
#define LBM_BND_PERIODIC_H
#include "bnd_cell.h"
#include "bnd_interface.h"
#include "lbm/constants.h"
#include "lbm/equilibrium_func.h"
#include "lbm/variables.h"

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBndCell_periodic : public LBMBndCell<LBTYPE> {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<LBEquationType::Navier_Stokes, NDIM>;

 public:
  LBMBndCell_periodic() = delete;
  LBMBndCell_periodic(const GInt cellId, const VectorD<dim(LBTYPE)>& _normal) : LBMBndCell<LBTYPE>(cellId, _normal) {}
  virtual ~LBMBndCell_periodic() = default;

  // deleted constructors not needed
  LBMBndCell_periodic(const LBMBndCell_periodic&)                    = default;
  LBMBndCell_periodic(LBMBndCell_periodic&&)                         = delete;
  auto operator=(const LBMBndCell_periodic&) -> LBMBndCell_periodic& = delete;
  auto operator=(LBMBndCell_periodic&&) -> LBMBndCell_periodic&      = delete;


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
        std::array<GDouble, NDIM> coord;
        for(GInt dir = 0; dir < NDIM; ++dir) {
          if(std::abs(normal()[dir]) > 0) {
            coord[dir] = bb.min(dir);
          } else {
            coord[dir] = center[dir] + LBMethod<LBTYPE>::m_dirs[dist][dir] * cellLength;
          }
        }
        // determine if this direction is within the bounding box and hence needs to be set
        if(bb.isInside(coord.data())) {
          m_bndIndex[m_noSetDists] = dist;
          ++m_noSetDists;
        }
      }
    }

    // todo: this needs to be moved since we cannot guarantee double matching cells are handled correctly
    for(GInt id = 0; id < m_noSetDists; ++id) {
      m_linkedCell[id] = INVALID_CELLID;
      const GInt dist  = m_bndIndex[id];

      // determine center of cell at one side of the periodic boundary condition
      auto centerA = center;
      for(GInt dir = 0; dir < dim(LBTYPE); ++dir) {
        if(std::abs(normal()[dir]) > 0) {
          continue;
        }
        centerA[dir] += LBMethod<LBTYPE>::m_dirs[dist][dir] * cellLength;
      }

      //      GDouble minDist = std::numeric_limits<GDouble>::max();
      static constexpr GDouble max_cell_match_dist = 10 * GDoubleEps;
      auto                     cellMatch           = [&](const auto& centerA_, const auto& centerB_) {
        for(GInt dir = 0; dir < dim(LBTYPE); ++dir) {
          const GDouble diff = std::abs(centerA_[dir] - centerB_[dir]);
          if(diff <= max_cell_match_dist) {
            return true;
          }
        }
        return false;
      };

      // connect cell of the boundary surface and connected surface
      for(const GInt cellIdB : surfConnected->getCellList()) {
        const auto& centerB = surfConnected->center(cellIdB);
        if(cellMatch(centerA, centerB)) {
          m_linkedCell[id] = cellIdB;
          break;
        }
      }

      debugMsg(surfConnected, centerA, id);

      ASSERT(m_linkedCell[id] != INVALID_CELLID,
             "No cell to link has been found for: " + std::to_string(mapped()) + " at: " + strStreamify<NDIM>(centerA).str());
    }
    errorCheck();
  }

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) {
    if(!std::isnan(m_pressure)) {
      // need to calculate a pressure drop across the periodic boundary
      for(GInt dist = 0; dist < NDIST; ++dist) {
        fold(m_linkedCell[0], dist) =
            eq::defaultEq<LBTYPE>(dist, m_pressure, &vars(mapped(), VAR::velocity(0))) + f(mapped(), dist) - feq(mapped(), dist);
      }
      vars(m_linkedCell[0], VAR::rho()) = m_pressure;
    } else {
      //    do the same as in the propagation step
      for(GInt id = 0; id < m_noSetDists; ++id) {
        const GInt dist              = m_bndIndex[id];
        fold(m_linkedCell[id], dist) = f(mapped(), dist);
      }
      vars(m_linkedCell[0], VAR::rho()) = 1.0;
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
             const std::function<GDouble&(GInt, GInt)>& /*vars*/) {
    // intentionally empty
  }


 private:
  using LBMBndCell<LBTYPE>::mapped;
  using LBMBndCell<LBTYPE>::normal;
  using LBMBndCell<LBTYPE>::init;

  void debugMsg(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surfConnected, const VectorD<NDIM>& centerA, const GInt id) {
    if(m_linkedCell[id] == INVALID_CELLID && DEBUG_LEVEL == Debug_Level::max_debug) {
      GDouble min       = std::numeric_limits<GDouble>::max();
      GInt    minCellId = INVALID_CELLID;
      ASSERT(!surfConnected->getCellList().empty(), "Surface is empty!");
      for(const GInt cellIdB : surfConnected->getCellList()) {
        const auto& centerB = surfConnected->center(cellIdB);
        for(GInt dir = 0; dir < dim(LBTYPE); ++dir) {
          const GDouble diff = std::abs(centerA[dir] - centerB[dir]);
          if(diff < min) {
            min       = diff;
            minCellId = cellIdB;
          } else if(minCellId == INVALID_CELLID) {
            cerr0 << cellIdB << std::endl;
          }
        }
      }
      cerr0 << "Only found a matching cell to: " << min << " with cellId: " << minCellId << std::endl;
      cerr0 << "For the location: " << strStreamify<NDIM>(centerA).str() << " with the matching cell center "
            << strStreamify<NDIM>(surfConnected->center(minCellId)).str() << std::endl;
      //      cerr0 << "The orientation was " << dist << std::endl;
      cerr0 << "CellLength " << surfConnected->cellLength(mapped()) << std::endl;
    }
  }

  void errorCheck() {
    // check if any cell has been linked twice
    for(GInt idA = 0; idA < m_noSetDists; ++idA) {
      for(GInt idB = idA + 1; idB < m_noSetDists; ++idB) {
        if(m_linkedCell[idA] == m_linkedCell[idB]) {
          TERMM(-1, "Invalid periodic bnd (cell has been linked twice)");
        }
      }
    }
  }

  std::array<GInt, noDists(LBTYPE)> m_bndIndex{};
  std::array<GInt, noDists(LBTYPE)> m_linkedCell{};
  GInt                              m_noSetDists = 0;
  GDouble                           m_pressure   = NAN;
};

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBnd_Periodic : public LBMBndInterface {
 public:
  LBMBnd_Periodic(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surfConnected,
                  const json& properties)
    : m_bnd(surf), m_bndConnected(surfConnected) {
    ASSERT(surfConnected->size() > 0, "Invalid connected surface");

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
  const SurfaceInterface*                               m_bnd          = nullptr;
  const SurfaceInterface*                               m_bndConnected = nullptr;
};
#endif // LBM_BND_PERIODIC_H
