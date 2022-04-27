#ifndef LPT_FORCES_H
#define LPT_FORCES_H

#include <sfcmm_common.h>

// Note: all forces are based on a=forces/mass so they are given as acceleration!!!!
namespace force {

/// Available force models
enum class Model {
  constDensityRatioGrav,              // all particles have the same density ratio to the ambient which is constant
  constDensityRatioGravStokesDrag,    // s.a. + Stokes Drag (Re_p << 1)
  constDensityRatioGravNlinDrag,      // s.a. + Non-linear Drag
  constDensityRatioGravBuo,           // s.a. but with non-neglible ambient density
  constDensityRatioGravBuoStokesDrag, // s.a. + Stokes Drag (Re_p << 1)
  constDensityRatioGravBuoNlinDrag,   // s.a. + Non-linear Drag
  constDensityaGravBuo,               // assume ambient density to be constant + Gravity and Buoancy
  constDensityaGravBuoStokesDrag,     // s.a. + + Stokes Drag (Re_p << 1)
  constDensityaGravBuoNlinDrag        // s.a. + Non-linear Drag
};

static constexpr std::array<std::string_view, 9> ModelName = {
    "constDensityRatioGrav",    "constDensityRatioGravStokesDrag",    "constDensityRatioGravNlinDrag",
    "constDensityRatioGravBuo", "constDensityRatioGravBuoStokesDrag", "constDensityRatioGravBuoNlinDrag",
    "constDensityaGravBuo",     "constDensityaGravBuoStokesDrag",     "constDensityaGravBuoNlinDrag"};

/// Convert between string and model type
/// \param modelName model type associated string
/// \return model type associated to the string
static constexpr auto model(const std::string_view modelName) -> Model {
  if(modelName == "constDensityRatioGrav") {
    return Model::constDensityRatioGrav;
  }
  if(modelName == "constDensityRatioGravStokesDrag") {
    return Model::constDensityRatioGravStokesDrag;
  }
  if(modelName == "constDensityRatioGravNlinDrag") {
    return Model::constDensityRatioGravNlinDrag;
  }
  if(modelName == "constDensityRatioGravBuo") {
    return Model::constDensityRatioGravBuo;
  }
  if(modelName == "constDensityRatioGravBuoStokesDrag") {
    return Model::constDensityRatioGravBuoStokesDrag;
  }
  if(modelName == "constDensityRatioGravBuoNlinDrag") {
    return Model::constDensityRatioGravBuoNlinDrag;
  }
  if(modelName == "constDensityaGravBuo") {
    return Model::constDensityaGravBuo;
  }
  if(modelName == "constDensityaGravBuoStokesDrag") {
    return Model::constDensityaGravBuoStokesDrag;
  }
  if(modelName == "constDensityaGravBuoNlinDrag") {
    return Model::constDensityaGravBuoNlinDrag;
  }

  TERMM(-1, "Invalid model configuration!");
}

/// Available drag models
enum class DragModel {
  Stokes,          // linear drag valid for Stokes flow (Rep << 1)
  SchillerNaumann, // non-linear drag by Schiller-Naumann
  PinskyKhain,     // non-linear drag by Pinsky&Khain, J. Aerosol Sci. 28(7), 1177-1214 (1997).
  Putnam61         // Mixture Formation in Internal Combustion Engines, 2005, Carsten Baumgarten
};

/// Calculate gravity
/// \tparam NDIM Dimensionality
/// \param gravity Gravity
/// \return Gravity
template <GInt NDIM>
static constexpr auto gravity(const VectorD<NDIM>& gravity) {
  return gravity;
}

/// Calculate net gravity i.e. gravity - buoyancy effect
/// \tparam NDIM Dimensionality
/// \param gravity Gravity
/// \param ambientDensity Ambient density of the particle surrounding medium
/// \param particleDensity Particle density
/// \return Net Gravity effect
template <GInt NDIM>
static constexpr auto gravityBuoyancy(const VectorD<NDIM>& gravity, const GDouble ambientDensity, const GDouble particleDensity) {
  return gravity * (1.0 - ambientDensity / particleDensity);
}

/// Calculate the drag force using the Stokes drag formula
/// \tparam NDIM Dimensionality
/// \param radius Radius of the particle
/// \param particleDensity Particle density
/// \param ambientViscosity Density of the particle surrounding medium
/// \param velocity Speed of the particle
/// \return Stokes drag force
template <GInt NDIM>
static constexpr auto dragForceStokes(const GDouble radius, const GDouble particleDensity, const GDouble ambientViscosity,
                                      const VectorD<NDIM>& velocity) {
  return 18.0 * ambientViscosity / (4.0 * radius * radius * particleDensity) * velocity;
}


/// Calculates the drag coefficient
/// \tparam DM Drag model to be used
/// \param rep Particle Reynoldsnumber
/// \return Drag coefficient
template <DragModel DM>
static constexpr auto dragCoefficient(const GDouble rep) {
  switch(DM) {
    case DragModel::Stokes:
      return 1.0;
    case DragModel::SchillerNaumann:
      return 1.0 + 0.15 * gcem::pow(rep, 0.687);
    case DragModel::PinskyKhain:
      return 1.0 + 0.17 * gcem::pow(rep, 0.632) + 1.0e-6 * gcem::pow(rep, 2.25);
    case DragModel::Putnam61:
      if(rep > 1000) {
        // Newton flow regime Re from 1k to 250k
        return 0.424 * rep / 24.0;
      } else if(rep <= 0.1) {
        // Stokes flow
        return 1.0;
      } else {
        // transition regime
        return 1.0 + gcem::pow(rep, 2.0 / 3.0) / 6.0;
      }
    default:
      TERMM(-1, "Invalid drag model!");
  }
}

} // namespace force

#endif // LPT_FORCES_H
