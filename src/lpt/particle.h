#ifndef LPT_PARTICLE_H
#define LPT_PARTICLE_H

#include <sfcmm_common.h>
#include "lpt/constants.h"

template <GInt NDIM, LPTType P>
class Particle {};

/// Base class for 1st order methods with one particle species
/// \tparam NDIM Dimensionality
template <GInt NDIM>
class Particle<NDIM, LPTType::Normal> {
 public:
  /// Number of double variables
  static constexpr std::byte m_noVars = static_cast<std::byte>(3 * NDIM + 5); // <= 112 bytes

  /// Index functions
  static constexpr auto center(const GInt dir) -> GInt { return dir; }
  static constexpr auto oldCenter(const GInt /*dir*/) -> GInt { return -1; }

  static constexpr auto velocity(const GInt dir) -> GInt { return NDIM + dir; }
  static constexpr auto oldVelocity(const GInt /*dir*/) -> GInt { return -1; }

  static constexpr auto a(const GInt dir) -> GInt { return 2 * NDIM + dir; }
  static constexpr auto oldA(const GInt /*dir*/) -> GInt { return -1; }

  static constexpr auto radius() -> GInt { return 3 * NDIM; }
  static constexpr auto oldRadius() -> GInt { return -1; }

  static constexpr auto density() -> GInt { return 3 * NDIM + 1; }

  static constexpr auto temperature() -> GInt { return 3 * NDIM + 2; }
  static constexpr auto oldTemperature() -> GInt { return -1; }

  static constexpr auto volume() -> GInt { return 3 * NDIM + 3; }
  static constexpr auto oldVolume() -> GInt { return -1; }

  static constexpr auto DC() -> GInt { return 3 * NDIM + 4; }
};

/// Base class for 2nd order and higher methods
/// \tparam NDIM Dimensionality
template <GInt NDIM>
class Particle<NDIM, LPTType::High> {
 public:
  /// Number of double variables
  static constexpr std::byte m_noVars = static_cast<std::byte>(6 * NDIM + 8); // <= 208 bytes

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

  static constexpr auto volume() -> GInt { return 6 * NDIM + 5; }
  static constexpr auto oldVolume() -> GInt { return 6 * NDIM + 6; }

  static constexpr auto DC() -> GInt { return 6 * NDIM + 7; }
};


/// Data access class to pass the result around
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <GInt NDIM, LPTType P>
class ParticleData : Particle<NDIM, P> {
 public:
  ParticleData(const GDouble* particleRef) : m_particleRef(particleRef) {}
  [[nodiscard]] inline constexpr auto center(const GInt dir) const -> GDouble { return m_particleRef[Particle<NDIM, P>::center(dir)]; }

  [[nodiscard]] inline constexpr auto velocity(const GInt dir) const -> GDouble { return m_particleRef[Particle<NDIM, P>::velocity(dir)]; }

  [[nodiscard]] inline constexpr auto a(const GInt dir) const -> GDouble { return m_particleRef[Particle<NDIM, P>::a(dir)]; }

  [[nodiscard]] inline constexpr auto radius() const -> GDouble { return m_particleRef[Particle<NDIM, P>::radius()]; }

  [[nodiscard]] inline constexpr auto density() const -> GDouble { return m_particleRef[Particle<NDIM, P>::density()]; }

  [[nodiscard]] inline constexpr auto temperature() const -> GDouble { return m_particleRef[Particle<NDIM, P>::temperature()]; }

  [[nodiscard]] inline constexpr auto volume() const -> GDouble { return m_particleRef[Particle<NDIM, P>::volume()]; }

  [[nodiscard]] inline constexpr auto DC() const -> GDouble { return m_particleRef[Particle<NDIM, P>::DC()]; }

 private:
  const GDouble* m_particleRef;
};

#endif // LPT_PARTICLE_H
