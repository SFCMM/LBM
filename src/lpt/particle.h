#ifndef LPT_PARTICLE_H
#define LPT_PARTICLE_H

#include <sfcmm_common.h>

enum class LPTType { Normal, High };

template <GInt NDIM, LPTType P>
class Particle {};

// Base class for 1st order methods with one particle species
template <GInt NDIM>
class Particle<NDIM, LPTType::Normal> {
 public:
  /// Number of double variables
  static constexpr std::byte m_noVars = static_cast<std::byte>(2 * NDIM + 3); // <= 72 bytes

  /// Index functions
  static constexpr auto coord(const GInt dir) -> GInt { return dir; }
  static constexpr auto oldCoord(const GInt dir) -> GInt { return -1; }

  static constexpr auto velocity(const GInt dir) -> GInt { return NDIM + dir; }
  static constexpr auto oldVelocity(const GInt dir) -> GInt { return -1; }

  static constexpr auto radius() -> GInt { return 2 * NDIM; }
  static constexpr auto oldRadius() -> GInt { return -1; }

  static constexpr auto density() -> GInt { return 2 * NDIM + 1; }

  static constexpr auto temperature() -> GInt { return 2 * NDIM + 2; }
  static constexpr auto oldTemperature() -> GInt { return -1; }
};

// Base class for 2nd order and higher methods
template <GInt NDIM>
class Particle<NDIM, LPTType::High> {
 public:
  /// Number of double variables
  static constexpr std::byte m_noVars = static_cast<std::byte>(4 * NDIM + 5); // <= 136 bytes

  /// Index functions
  static constexpr auto coord(const GInt dir) -> GInt { return dir; }
  static constexpr auto oldCoord(const GInt dir) -> GInt { return NDIM + dir; }

  static constexpr auto velocity(const GInt dir) -> GInt { return 2 * NDIM + dir; }
  static constexpr auto oldVelocity(const GInt dir) -> GInt { return 3 * NDIM + dir; }

  static constexpr auto radius() -> GInt { return 4 * NDIM; }
  static constexpr auto oldRadius() -> GInt { return 4 * NDIM + 1; }

  static constexpr auto density() -> GInt { return 2 * NDIM + 2; }

  static constexpr auto temperature() -> GInt { return 2 * NDIM + 3; }
  static constexpr auto oldTemperature() -> GInt { return 2 * NDIM + 4; }
};

#endif // LPT_PARTICLE_H
