#ifndef LPT_CONSTANTS_H
#define LPT_CONSTANTS_H

#include "common/term.h"

// default values
static constexpr GInt default_number_particles_capacity = 100;

// Type which switches between first and second order accuracy
enum class LPTType { Normal, High };

/// Initial conditions
// todo: implement
enum class LPTInitCond { none, load_from_CSV, randomvol_pos, randomplane_pos };

/// Initial condition names as string
static constexpr std::array<std::string_view, 4> LPTInitCondName = {"none", "load_from_CSV", "randomvol_pos", "randomplane_pos"};

/// Convert between string an LPTInitCond type
/// \param condName String name of the LPT Initial condition
/// \return LPTInitCond Type associated with the string name
static constexpr auto initCond(const std::string_view condName) -> LPTInitCond {
  if(condName == "none") {
    return LPTInitCond::none;
  }
  if(condName == "load_from_CSV") {
    return LPTInitCond::load_from_CSV;
  }
  if(condName == "randomvol_pos") {
    return LPTInitCond::randomvol_pos;
  }
  if(condName == "randomplane_pos") {
    return LPTInitCond::randomplane_pos;
  }
  TERMM(-1, "Invalid initial condition configuration!");
}

/// Available integration methods
// todo: implement
enum class IntegrationMethod {
  ForwardEuler,
  ForwardEulerPredCor, // s.a. + Predictor-Corrector Step
  ImplicitEuler,       // for high relative velocities (Note: underpredicts actual velocity)
  Berger20, // Implicit First-order Exponential Integrator Method as in Sven Berger et al., Large-Eddy Simulation Study of Biofuel Injection
            // in an Optical Direct Injection Engine (Note: stable and more accurate for high-relative velocities)
  Berger20PredCor // s.a. + Predictor-Corrector Step
};

/// Method that generate particles
// todo:implement
enum class GenerationMethod { None, ConstantRate, InjectionModel };

/// String names of generationMethod
static constexpr std::array<std::string_view, 4> GenerationMethodName = {"none", "constant_rate", "injection"};

/// Convert between string and GenerationMethod type
/// \param generationName GenerationMethod type associated string
/// \return GenerationMethod type associated to the string
static constexpr auto generationMethod(const std::string_view generationName) -> GenerationMethod {
  if(generationName == "none") {
    return GenerationMethod::None;
  }
  if(generationName == "constant_rate") {
    return GenerationMethod::ConstantRate;
  }
  if(generationName == "injection") {
    return GenerationMethod::InjectionModel;
  }

  TERMM(-1, "Invalid generation method configuration!");
}
#endif // LPT_CONSTANTS_H
