#ifndef LPT_SOLVER_H
#define LPT_SOLVER_H

#include <Eigen/Core>
#include <sfcmm_common.h>
#include "cartesiangrid.h"
#include "common/configuration.h"
#include "common/random.h"
#include "interface/solver_interface.h"
#include "particle.h"

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
  void output(const GBool forced);

  inline auto velocity(const GInt pid) -> Eigen::Map<VectorD<NDIM>> {
    return Eigen::Map<VectorD<NDIM>>(&m_vars[pid * NVARS + PTYPE::velocity(0)]);
  }

  inline auto center(const GInt pid) -> Eigen::Map<VectorD<NDIM>> {
    return Eigen::Map<VectorD<NDIM>>(&m_vars[pid * NVARS + PTYPE::center(0)]);
  }

  inline auto a(const GInt pid) -> Eigen::Map<VectorD<NDIM>> { return Eigen::Map<VectorD<NDIM>>(&m_vars[pid * NVARS + PTYPE::a(0)]); }

  inline auto density(const GInt pid) -> GDouble& { return m_vars[pid * NVARS + PTYPE::density()]; }

  inline auto volume(const GInt pid) -> GDouble& { return m_vars[pid * NVARS + PTYPE::volume()]; }

  inline auto radius(const GInt pid) -> GDouble& { return m_vars[pid * NVARS + PTYPE::radius()]; }

  inline auto temperature(const GInt pid) -> GDouble& { return m_vars[pid * NVARS + PTYPE::temperature()]; }


  /// Configuration
  GString     m_exe;
  GString     m_configurationFileName;
  GBool       m_benchmark        = false;
  GInt        m_capacity         = 1000;
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
  GInt                 m_timeStep = 0;
  GDouble              m_dt       = NAN;
};


#endif // LPT_SOLVER_H
