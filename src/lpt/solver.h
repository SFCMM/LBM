#ifndef LPT_SOLVER_H
#define LPT_SOLVER_H

#include <Eigen/Core>
#include <sfcmm_common.h>
#include "cartesiangrid.h"
#include "common/configuration.h"
#include "common/random.h"
#include "forces.h"
#include "interface/solver_interface.h"
#include "particle.h"

enum class IntegrationMethod {
  ForwardEuler,
  ForwardEulerPredCor, // s.a. + Predictor-Corrector Step
  ImplicitEuler,       // for high relative velocities (Note: underpredicts actual velocity)
  Berger20, // Implicit First-order Exponential Integrator Method as in Sven Berger et al., Large-Eddy Simulation Study of Biofuel Injection
            // in an Optical Direct Injection Engine (Note: stable and more accurate for high-relative velocities)
  Berger20PredCor // s.a. + Predictor-Corrector Step
};

/// Solver for Lagrange particles
/// Memory consumption 1MB = ~10000 Particles (Normal), ~5000 Particles (Evaporation, High Accuracy)
/// \tparam DEBUG_LEVEL
/// \tparam NDIM
/// \tparam P
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
class LPTSolver : public Runnable, private Configuration, private RandomGenerator {
 private:
  // Type declarations and template variables only
  using PTYPE                 = Particle<NDIM, P>;
  static constexpr GInt NVARS = static_cast<GInt>(PTYPE::m_noVars);

 public:
  LPTSolver(GInt32 domainId, GInt32 noDomains) : m_domainId(domainId), m_noDomains(noDomains){};
  ~LPTSolver() override       = default;
  LPTSolver(const LPTSolver&) = delete;
  LPTSolver(LPTSolver&&)      = delete;
  auto operator=(const LPTSolver&) -> LPTSolver& = delete;
  auto operator=(LPTSolver&&) -> LPTSolver& = delete;

  void               init(int argc, GChar** argv, GString config_file) override;
  void               initBenchmark(int argc, GChar** argv) override;
  auto               run() -> GInt override;
  [[nodiscard]] auto grid() const -> const CartesianGrid<DEBUG_LEVEL, NDIM>& override { TERMM(-1, "Not implemented"); };
  void               transferGrid(const GridInterface& /*grid*/) override { TERMM(-1, "Not implemented"); };

  /// Memory used in Kbytes
  /// \return KBytes used
  [[nodiscard]] auto mem() const -> GDouble { return m_vars.size() * 8.0 / DKBIT; }

 private:
  void init(int argc, GChar** argv);
  void initTimers();
  void allocateMemory();
  void loadConfiguration();

  void initialCondition();
  void init_randomVolPos();

  void timeStep();
  template <force::Model FM, IntegrationMethod IM>
  void calcA();
  template <IntegrationMethod IM>
  void timeIntegration();

  void compareToAnalyticalResult();
  void output(const GBool forced, const GString& postfix = "");

  // Variable accessor functions
  inline auto part(const GInt pid) -> ParticleData<NDIM, P> { return ParticleData<NDIM, P>(&m_vars[pid * NVARS]); }

  /// Set/Get velocity of a particle
  /// \param pid Particle Id
  /// \return Velocity of the particle
  inline auto velocity(const GInt pid) -> Eigen::Map<VectorD<NDIM>> {
    return Eigen::Map<VectorD<NDIM>>(&m_vars[pid * NVARS + PTYPE::velocity(0)]);
  }

  /// Set/Get velocity of a particle
  /// \param pid Particle Id
  /// \param dir Direction
  /// \return Velocity of the particle
  inline auto velocity(const GInt pid, const GInt dir) -> GDouble& { return m_vars[pid * NVARS + PTYPE::velocity(dir)]; }

  /// Set/Get center of a particle
  /// \param pid Particle Id
  /// \return Center of the particle
  inline auto center(const GInt pid) -> Eigen::Map<VectorD<NDIM>> {
    return Eigen::Map<VectorD<NDIM>>(&m_vars[pid * NVARS + PTYPE::center(0)]);
  }

  /// Set/Get acceleration of a particle
  /// \param pid Particle Id
  /// \return acceleration of the particle
  inline auto a(const GInt pid) -> Eigen::Map<VectorD<NDIM>> { return Eigen::Map<VectorD<NDIM>>(&m_vars[pid * NVARS + PTYPE::a(0)]); }

  /// Set/Get density of a particle
  /// \param pid Particle Id
  /// \return density of the particle
  inline auto density(const GInt pid) -> GDouble& { return m_vars[pid * NVARS + PTYPE::density()]; }

  /// Set/Get volume of a particle
  /// \param pid Particle Id
  /// \return volume of the particle
  inline auto volume(const GInt pid) -> GDouble& { return m_vars[pid * NVARS + PTYPE::volume()]; }

  /// Set/Get radius of a particle
  /// \param pid Particle Id
  /// \return radius of the particle
  inline auto radius(const GInt pid) -> GDouble& { return m_vars[pid * NVARS + PTYPE::radius()]; }

  /// Set/Get temperature of a particle
  /// \param pid Particle Id
  /// \return temperature of the particle
  inline auto temperature(const GInt pid) -> GDouble& { return m_vars[pid * NVARS + PTYPE::temperature()]; }

  /// Set/Get drag coefficient of a particle
  /// \param pid Particle Id
  /// \return drag coefficient of the particle
  inline auto DC(const GInt pid) -> GDouble& { return m_vars[pid * NVARS + PTYPE::DC()]; }


  /// Configuration
  GString     m_exe;
  GString     m_configurationFileName;
  GBool       m_benchmark        = false;
  GInt        m_capacity         = default_number_particles_capacity;
  LPTInitCond m_initialCondition = LPTInitCond::none;

  /// MPI
  GInt32 m_domainId  = -1;
  GInt32 m_noDomains = -1;

  /// Output
  GString m_outputDir              = "out/";
  GString m_solutionFileName       = "solution";
  GInt    m_outputInfoInterval     = 100;
  GInt    m_outputSolutionInterval = 1000;
  GBool   m_generatePath           = true;

  /// Variables
  std::vector<GDouble> m_vars;
  GInt                 m_timeStep    = 0;
  GDouble              m_currentTime = 0;
  GInt                 m_noParticles = 0;
  GDouble              m_dt          = NAN;

  VectorD<NDIM> m_init_v;

  // ambient condition
  //   VectorD<NDIM>        m_gravity     = NAN_LIST<NDIM>();
  VectorD<NDIM> m_gravity;
  GDouble       m_rho_a_infty = NAN;
  GDouble       m_nu_a_infty  = NAN;
  //  VectorD<NDIM> m_velo_a_infty =  NAN_LIST<NDIM>();
  VectorD<NDIM> m_velo_a_infty;
};


#endif // LPT_SOLVER_H
