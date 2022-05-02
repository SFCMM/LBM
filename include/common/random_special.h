// SPDX-License-Identifier: BSD-3-Clause

#ifndef COMMON_RANDOM_SPECIAL_H
#define COMMON_RANDOM_SPECIAL_H

#include "common/math/mathfunctions.h"
#include "gcem.hpp"
#include "geometry/circle.h"
#include "randxor.h"
#include "sfcmm_types.h"

enum class Random_Dist { none, normal };

/// Generate a random normalized direction within an cone defined by its center axis orientation and an opening angle.
/// \param centerAxis Center orientation of the cone
/// \param openingAngle Opening angle of the cone in degrees
/// \param dist Distribution used (default: uniform)
/// \param distAttrb Distribution attribute (default: 0)
/// \return A random normalized direction within the defined cone.
inline auto randomNormalizedDirection_inCone(const VectorD<3>& coneAxis, const GDouble openingAngle, randxor& prng,
                                             const Random_Dist /*dist*/ = Random_Dist::none, const GDouble /*distAttrb*/ = 0.0)
    -> VectorD<3> {
  // angle between center cone axis and cone shell
  GDouble distedAngle = 0;
  // todo: implement
  //  if(dist == Random_Dist::normal) {
  //    //~68% of all particles are within in the center cone for distCoeff = 1.0
  //    normal_distribution<GDouble> distAngle(0, distAttrb * openingAngle);
  //    // reject angles that are too large
  //    do {
  //      distedAngle = distAngle(PRNG);
  //    } while(distedAngle > openingAngle);
  //  } else {
  // uniform
  distedAngle = openingAngle;
  //  }
  const GDouble coneAngleRad = toRadians(distedAngle);

  // 1. normalize coneAxis
  VectorD<3> centerAxis = coneAxis.normalized();

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

  if(std::abs(coneAxis[0]) < std::numeric_limits<GDouble>::epsilon()) {
    v1[0] = 1.0;
  } else if(std::abs(coneAxis[1]) < std::numeric_limits<GDouble>::epsilon()) {
    v1[1] = 1.0;
  } else if(std::abs(coneAxis[2]) < std::numeric_limits<GDouble>::epsilon()) {
    v1[2] = 1.0;
  }

  // 2.2 Use this vector to determine the first basis vector by calculating v1 = v1 x centerAxis
  v1 = v1.cross(centerAxis);

  // 2.3 Calculate crossproduct of v1 and centerAxis to determine the last vector v2
  VectorD<3> v2 = v1.cross(centerAxis);


  // 2.4 Normalize v1 and v2
  v1 = v1.normalized();
  v2 = v2.normalized();

  // random point on unit circle
  const GDouble phi = prng.double_value(0, 2 * M_PI);

  // random angle within the given opening angle
  const GDouble omega = gcem::acos(prng.double_value(gcem::cos(coneAngleRad), 1.0));

  VectorD<3> tmp = gcem::sin(omega) * (gcem::cos(phi) * v1 + gcem::sin(phi) * v2) + gcem::cos(omega) * centerAxis;

  // 3. Use orthonormal basis to determine points on circle
  return tmp.normalized();
}

/// Obtain a random point in a circle on the given plane by the normal around the center (0.0, 0.0, 0.0).
/// \param[in] normal The plane's normal vector.
/// \param[in] diameter Diameter of the circle.
/// \param[in] prng Random number generator reference
/// \return Point on the defined circle.
inline auto randomPoint_inCircle(const VectorD<3>& normal, const GDouble diameter, randxor& prng) -> VectorD<3> {
  const GDouble random_angle    = prng.double_value(0, 2 * M_PI);
  const GDouble random_diameter = prng.double_value() * diameter;
  return pointOnCircleShell(normal, random_diameter, random_angle);
}

/// Obtain a random point on a circular shell with the given plane by the normal  around the center (0.0, 0.0, 0.0).
/// \param[in] normal The plane's normal vector.
/// \param[in] diameter Diameter of the circle.
/// \param[in] prng Random number generator reference
/// \return Point on the defined circle.
inline auto randomPoint_onCircleShell(const VectorD<3>& normal, const GDouble diameter, randxor& prng) -> VectorD<3> {
  const GDouble random_angle = prng.double_value(0, 2 * M_PI + GDoubleEps);
  return pointOnCircleShell(normal, diameter, random_angle);
}

/// Generate a random value using the Rosin-Rammler distribution. This often used for the size distribution of droplets during injection.
/// \tparam RR_SPREAD Constant spread factor
/// \param mean Mean value of the distribution
/// \param prng Random number generator reference to be used.
/// \return Random value distributed by the Rosin-Rammler distribution within [0.3 * mean, 3.0 * mean]
template <GInt RR_SPREAD = 3>
inline constexpr auto randomDist_rosinRammler(const GDouble mean, randxor& prng) -> GDouble {
  GDouble rr_max = 3.0 * mean;
  GDouble rr_min = 0.3 * mean;

  const GDouble K = 1.0 - gcem::exp(-gcem::pow((rr_max - rr_min) / mean, RR_SPREAD));
  const GDouble x = prng.double_value();

  return rr_min + mean * gcem::pow(-gcem::log(1.0 - x * K), 1.0 / RR_SPREAD);
}

#endif // COMMON_RANDOM_SPECIAL_H
