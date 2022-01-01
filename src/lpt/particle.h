#ifndef LPT_PARTICLE_H
#define LPT_PARTICLE_H

#include <sfcmm_common.h>
#include "lpt/constants.h"


template <GInt NDIM, LPTType P>
class Particle {};

// Base class for 1st order methods with one particle species
template <GInt NDIM>
class Particle<NDIM, LPTType::Normal> {
 public:
  /// Number of double variables
  static constexpr std::byte m_noVars = static_cast<std::byte>(3 * NDIM + 3); // <= 96 bytes

  /// Index functions
  static constexpr auto center(const GInt dir) -> GInt { return dir; }
  static constexpr auto oldCenter(const GInt dir) -> GInt { return -1; }

  static constexpr auto velocity(const GInt dir) -> GInt { return NDIM + dir; }
  static constexpr auto oldVelocity(const GInt dir) -> GInt { return -1; }

  static constexpr auto a(const GInt dir) -> GInt { return 2 * NDIM + dir; }
  static constexpr auto oldA(const GInt dir) -> GInt { return -1; }

  static constexpr auto radius() -> GInt { return 3 * NDIM; }
  static constexpr auto oldRadius() -> GInt { return -1; }

  static constexpr auto density() -> GInt { return 3 * NDIM + 1; }

  static constexpr auto temperature() -> GInt { return 3 * NDIM + 2; }
  static constexpr auto oldTemperature() -> GInt { return -1; }
};

// Base class for 2nd order and higher methods
template <GInt NDIM>
class Particle<NDIM, LPTType::High> {
 public:
  /// Number of double variables
  static constexpr std::byte m_noVars = static_cast<std::byte>(6 * NDIM + 5); // <= 184 bytes

  /// Index functions
  static constexpr auto coord(const GInt dir) -> GInt { return dir; }
  static constexpr auto oldCoord(const GInt dir) -> GInt { return NDIM + dir; }

  static constexpr auto velocity(const GInt dir) -> GInt { return 2 * NDIM + dir; }
  static constexpr auto oldVelocity(const GInt dir) -> GInt { return 3 * NDIM + dir; }

  static constexpr auto a(const GInt dir) -> GInt { return 4 * NDIM + dir; }
  static constexpr auto oldA(const GInt dir) -> GInt { return 5 * NDIM + dir; }

  static constexpr auto radius() -> GInt { return 6 * NDIM; }
  static constexpr auto oldRadius() -> GInt { return 6 * NDIM + 1; }

  static constexpr auto density() -> GInt { return 6 * NDIM + 2; }

  static constexpr auto temperature() -> GInt { return 6 * NDIM + 3; }
  static constexpr auto oldTemperature() -> GInt { return 6 * NDIM + 4; }
};

#endif // LPT_PARTICLE_H
