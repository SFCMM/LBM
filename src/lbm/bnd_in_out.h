#ifndef LBM_BND_IN_OUT_H
#define LBM_BND_IN_OUT_H
#include "bnd_interface.h"
#include "common/surface.h"
#include "constants.h"
#include "equilibrium_func.h"
#include "variables.h"

// Antibounceback pressure boundary condition
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBnd_Pressure : public LBMBndInterface {
 private:
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;

  //  static constexpr GInt NVAR = noVars<LBTYPE>(LBEquationType::Navier_Stokes);

  using VAR = LBMVariables<LBEquationType::Navier_Stokes, NDIM>;

 public:
  LBMBnd_Pressure(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties)
    : m_bnd(surf), m_pressure(config::required_config_value<GDouble>(properties, "pressure")) {}
  ~LBMBnd_Pressure() override = default;

  // deleted constructors not needed
  LBMBnd_Pressure(const LBMBnd_Pressure&) = delete;
  LBMBnd_Pressure(LBMBnd_Pressure&&)      = delete;
  auto operator=(const LBMBnd_Pressure&) -> LBMBnd_Pressure& = delete;
  auto operator=(LBMBnd_Pressure&&) -> LBMBnd_Pressure& = delete;

  void initCnd(const std::function<GDouble&(GInt, GInt)>& vars) override {
    for(const auto bndCellId : m_bnd->getCellList()) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(bndCellId, VAR::rho()) = m_pressure;
      }
    }
  }

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& vars) override {
    // value is set first so that this value can be used in other boundary conditions
    // todo: this is wrong if there are two pressure boundary conditions with different values
    for(const auto bndCellId : m_bnd->getCellList()) {
      vars(bndCellId, VAR::rho()) = m_pressure;
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>&    fold,
             const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& vars) override {
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

    for(const auto bndCellId : m_bnd->getCellList()) {
      const GDouble* normal          = m_bnd->normal_p(bndCellId);
      const GInt     insideDir       = extrapolationDir(normal);
      const GInt     insideNeighbor  = m_bnd->neighbor(bndCellId, insideDir);
      const GInt     insideNeighbor2 = m_bnd->neighbor(insideNeighbor, insideDir);

      //      cerr0 << "insideDir " << insideDir << " " << strStreamify<NDIM>(VectorD<NDIM>(normal)).str() << " pressure " << m_pressure
      //            << std::endl;
      const VectorD<NDIM> V1(&vars(insideNeighbor, VAR::velocity(0)));
      const VectorD<NDIM> V2(&vars(insideNeighbor2, VAR::velocity(0)));
      //      cerr0 << "V1 before " << strStreamify<NDIM>(V1).str() << std::endl;
      //      cerr0 << "insideNGhbriD " << insideNeighbor << " of " << bndCellId << std::endl;
      ASSERT(insideNeighbor != INVALID_CELLID, "Invalid cell accessed");
      ASSERT(insideNeighbor2 != INVALID_CELLID, "Invalid cell accessed");

      VectorD<NDIM> extrpV = 1.5 * V1 - 0.5 * V2;
      //      VectorD<NDIM> extrpV = 0.5 * V1 + 0.5 * V2;
      //      cerr0 << "wallv after " << strStreamify<NDIM>(extrpV).str() << std::endl;

      //      for(GInt dir = 0; dir< NDIM; ++dir){
      //        vars(bndCellId, VAR::velocity(dir)) = extrpV[dir];
      //      }

      for(GInt dist = 0; dist < NDIST; ++dist) {
        // todo: actually depends on the type of bnd. this is only correct for walls!
        //  on a corner
        //        if(m_bnd->neighbor(bndCellId, 3) == INVALID_CELLID || m_bnd->neighbor(bndCellId, 2) == INVALID_CELLID) {
        //          extrpV.fill(0);
        //        }


        //        // recalculate equilibrium distribution
        //        cerr0 << "fold before " << fold(bndCellId, dist) << std::endl;
        //        fold(bndCellId, dist) = eq::defaultEq<LBTYPE>(dist, m_pressure, extrpV.data());
        //        //          fold(bndCellId, oppositeDist) -= 2 * eq::defaultEq<LBTYPE>(dist, m_pressure, extrpV.data());//totally wrong
        //        cerr0 << "fold after " << fold(bndCellId, dist) << std::endl;
        //
        //        vars(bndCellId, VAR::rho()) = m_pressure;
        //        vars(bndCellId, VAR::velocity(0)) = extrpV[0];
        //        vars(bndCellId, VAR::velocity(1)) = extrpV[1];

        // skip setting dists that are normal to the boundary
        if(dist < NDIST - 1 && m_bnd->neighbor(bndCellId, dist) == INVALID_CELLID
           && inDirection<NDIM>(VectorD<NDIM>(m_bnd->normal_p(bndCellId)), LBMethod<LBTYPE>::m_dirs[dist])) {
          // dist in reflected/opposite direction i.e. to the inside of the bnd
          GInt oppositeDist = LBMethod<LBTYPE>::oppositeDist(dist);
          //          cerr0 << "setting " << oppositeDist << " with " << dist << std::endl;

          // anti bounceback i.e. distribution that hits the wall is reflected to the opposite distribution direction with a negative sign
          fold(bndCellId, oppositeDist) = -fpre(bndCellId, dist) + 2 * eq::symmEq<LBTYPE>(dist, m_pressure, extrpV.data());

          //          cerr0 << "fold " << fold(bndCellId, oppositeDist) << " fpre " << fpre(bndCellId, dist) << std::endl;

          //          const GDouble currentPressure = vars(bndCellId, VAR::rho());
          //          const GDouble dP              = m_pressure - currentPressure;
          //          if(std::abs(dP) > m_pressure) {
          //            TERMM(-1, "testing");
          //          }
          //          cerr0 << "dp " << dP << " current " << currentPressure << std::endl;
        }
      }
    }
  }


 private:
  const SurfaceInterface* m_bnd = nullptr;

  GDouble m_pressure = 1.0;
};

#endif // LBM_BND_IN_OUT_H
