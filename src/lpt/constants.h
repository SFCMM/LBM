#ifndef LPT_CONSTANTS_H
#define LPT_CONSTANTS_H

// default values
static constexpr GInt default_number_particles_capacity = 100;

enum class LPTType { Normal, High };

enum class LPTInitCond { none, load_from_CSV, randomvol_pos, randomplane_pos };

enum class IntegrationMethod {
  ForwardEuler,
  ForwardEulerPredCor, // s.a. + Predictor-Corrector Step
  ImplicitEuler,       // for high relative velocities (Note: underpredicts actual velocity)
  Berger20, // Implicit First-order Exponential Integrator Method as in Sven Berger et al., Large-Eddy Simulation Study of Biofuel Injection
            // in an Optical Direct Injection Engine (Note: stable and more accurate for high-relative velocities)
  Berger20PredCor // s.a. + Predictor-Corrector Step
};

enum class GenerationMethod { None, ConstantRate, InjectionModel };

static constexpr std::array<std::string_view, 4> LPTInitCondName = {"none", "load_from_CSV", "randomvol_pos", "randomplane_pos"};

static constexpr auto getLPTInitCond(const std::string_view condName) -> LPTInitCond {
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
#endif // LPT_CONSTANTS_H
