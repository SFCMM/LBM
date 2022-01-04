#ifndef LPT_FORCES_H
#define LPT_FORCES_H

#include <sfcmm_common.h>

// Note: all forces are based on a=forces/mass so they are given as acceleration!!!!
namespace force {

enum class Model {
  constDensityRatioGrav,              // all particles have the same density ratio to the ambient which is constant
  constDensityRatioGravStokesDrag,    // s.a. + Stokes Drag (Re_p << 1)
  constDensityRatioGravBuo,           // s.a. but with non-neglible ambient density
  constDensityRatioGravBuoStokesDrag, // s.a. + Stokes Drag (Re_p << 1)
  constDensityaGravBuo,               // assume ambient density to be constant + Gravity and Buoancy
  constDensityaGravBuoStokesDrag      // s.a. + + Stokes Drag (Re_p << 1)
};

template <GInt NDIM>
static constexpr auto gravity(const VectorD<NDIM>& gravity) {
  return gravity;
}

template <GInt NDIM>
static constexpr auto gravityBuoyancy(const VectorD<NDIM>& gravity, const GDouble ambientDensity, const GDouble particleDensity) {
  return gravity * (1.0 - ambientDensity / particleDensity);
}

template <GInt NDIM>
static constexpr auto dragForceStokes(const GDouble radius, const GDouble particleDensity, const GDouble ambientViscosity,
                                      const VectorD<NDIM>& velocity) {
  return 18.0 * ambientViscosity / (4.0 * radius * radius * particleDensity) * velocity;
}

} // namespace force

#endif // LPT_FORCES_H
