// SPDX-License-Identifier: BSD-3-Clause

#ifndef COMMON_CIRCLE_H
#define COMMON_CIRCLE_H
#include <Eigen/Geometry>
#include "gcem.hpp"
#include "../sfcmm_types.h"

/// Obtain a point on a circular shell with the given plane around the center (0.0, 0.0, 0.0) given a diameter and angle.
/// Note: Rotation is anticlockwise but axes order is not guaranteed
/// \param[in] normal The plane's normal vector.
/// \param[in] diameter Diameter of the circle.
/// \param[in] phi Angle of the point
/// \return Point on the defined circle.
inline auto pointOnCircleShell(const VectorD<3>& normal, const GDouble diameter, const GDouble phi) -> VectorD<3> {
  // 1. normalize centerAxis
  VectorD<3> centerAxis = normal.normalized();

  // 2. find orthonormal basis to centerAxis
  // 2.1 choose a linear independent vector to centerAxis that lies in the plane a v1 * centerAxis = 1
  VectorD<3> v1({0.0, 0.0, 0.0});

  if(std::abs(centerAxis[0]) > 0.0) {
    v1[0] = 1.0 / centerAxis[0];
  } else if(std::abs(centerAxis[1]) > 0.0) {
    v1[1] = 1.0 / centerAxis[1];
  } else {
    v1[2] = 1.0 / centerAxis[2];
  }

  if(std::abs(normal[0]) < std::numeric_limits<GDouble>::epsilon()) {
    v1[0] = 1.0;
  } else if(std::abs(normal[1]) < std::numeric_limits<GDouble>::epsilon()) {
    v1[1] = 1.0;
  } else if(std::abs(normal[2]) < std::numeric_limits<GDouble>::epsilon()) {
    v1[2] = 1.0;
  }

  // 2.2 Use this vector to determine the first basis vector by calculating v1 = v1 x centerAxis
  v1 = v1.cross(centerAxis);

  // 2.3 Calculate crossproduct of v1 and centerAxis to determine the last vector v2
  VectorD<3> v2 = v1.cross(centerAxis);


  // 2.4 Normalize v1 and v2
  v1 = v1.normalized();
  v2 = v2.normalized();

  // 3. Use orthonormal basis to determine points on circle
  return 0.5 * diameter * (gcem::cos(phi) * v1 + gcem::sin(phi) * v2);
}

#endif // COMMON_CIRCLE_H
