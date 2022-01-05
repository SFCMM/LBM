#ifndef LBM_SPHERE_H
#define LBM_SPHERE_H

#include <sfcmm_common.h>

namespace sphere {

/// Area of a sphere given a radius
/// \param radius Radius of the sphere
/// \return Area of the sphere
static constexpr auto areaR(const GDouble radius) -> GDouble { return 4.0 * PI * gcem::pow(radius, 2); }

/// Volume of a sphere given a radius
/// \param radius Radius of the sphere
/// \return Volume of the sphere
static constexpr auto volumeR(const GDouble radius) -> GDouble { return 4.0 / 3.0 * PI * gcem::pow(radius, 3); }

} // namespace sphere

#endif // LBM_SPHERE_H
