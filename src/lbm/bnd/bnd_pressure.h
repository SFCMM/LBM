#ifndef LBM_BND_PRESSURE_H
#define LBM_BND_PRESSURE_H
#include "bnd_interface.h"
#include "common/surface.h"
#include "lbm/constants.h"
#include "lbm/equilibrium_func.h"
#include "lbm/variables.h"

/// Antibounceback pressure boundary condition
/// Assumes currently that the density and temperature are constant which means pressure = density
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
  LBMBnd_Pressure(const LBMBnd_Pressure&)                    = delete;
  LBMBnd_Pressure(LBMBnd_Pressure&&)                         = delete;
  auto operator=(const LBMBnd_Pressure&) -> LBMBnd_Pressure& = delete;
  auto operator=(LBMBnd_Pressure&&) -> LBMBnd_Pressure&      = delete;


  void initCnd(const std::function<GDouble&(GInt, GInt)>& vars) override {
    for(const auto bndCellId : m_bnd->getCellList()) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(bndCellId, VAR::rho()) = m_pressure;
      }
    }
  }

  /// Set the pressure value (density) in the boundary cell
  /// \param vars
  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& vars) override {
    // value is set first so that this value can be used in other boundary conditions
    // todo: this is wrong if there are two pressure boundary conditions with different values
    for(const auto bndCellId : m_bnd->getCellList()) {
      vars(bndCellId, VAR::rho()) = m_pressure;
    }
  }

  /// Set the boundary condition
  /// \param fpre
  /// \param fold
  /// \param vars
  void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>&    fold,
             const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& vars) override {
    // Determine the inside direction
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
      const GDouble* normal    = m_bnd->normal_p(bndCellId);
      const GInt     insideDir = extrapolationDir(normal);
      // direct neighbor in the inside from the boundary cell
      const GInt directInsideNghbr = m_bnd->neighbor(bndCellId, insideDir);
      // inside neighbor of the boundary cell 2 cells to the inside
      const GInt insideNormalNghbr2 = m_bnd->neighbor(directInsideNghbr, insideDir);

      const VectorD<NDIM> directInsideNghbr_velocity(&vars(directInsideNghbr, VAR::velocity(0)));
      const VectorD<NDIM> insideNormalNghbr2_velocity(&vars(insideNormalNghbr2, VAR::velocity(0)));

      VectorD<NDIM> extrpV = 1.5 * directInsideNghbr_velocity - 0.5 * insideNormalNghbr2_velocity;

      for(GInt dist = 0; dist < NDIST; ++dist) {
        // todo: actually depends on the type of bnd. this is only correct for walls!
        //  on a corner
        //        if(m_bnd->neighbor(bndCellId, 3) == INVALID_CELLID || m_bnd->neighbor(bndCellId, 2) == INVALID_CELLID) {
        //          extrpV.fill(0);
        //        }


        vars(bndCellId, VAR::rho())       = m_pressure;
        vars(bndCellId, VAR::velocity(0)) = extrpV[0];
        vars(bndCellId, VAR::velocity(1)) = extrpV[1];

        // skip setting dists that are normal to the boundary
        if(dist < NDIST - 1 && m_bnd->neighbor(bndCellId, dist) == INVALID_CELLID
           && inDirection<NDIM>(VectorD<NDIM>(m_bnd->normal_p(bndCellId)), LBMethod<LBTYPE>::m_dirs[dist])) {
          // dist in reflected/opposite direction i.e. to the inside of the bnd
          GInt oppositeDist = LBMethod<LBTYPE>::oppositeDist(dist);

          // anti bounceback i.e. distribution that hits the wall is reflected to the opposite distribution direction with a negative sign
          fold(bndCellId, oppositeDist) = -fpre(bndCellId, dist) + 2 * eq::symmEq<LBTYPE>(dist, m_pressure, extrpV.data());
        }
      }
    }
  }


 private:
  const SurfaceInterface* m_bnd = nullptr;

  GDouble m_pressure = 1.0;
};

#endif // LBM_BND_PRESSURE_H
