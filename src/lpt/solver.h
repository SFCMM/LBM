#ifndef LPT_SOLVER_H
#define LPT_SOLVER_H

#include <sfcmm_common.h>
#include "common/configuration.h"
#include "interface/solver_interface.h"
#include "particle.h"

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
class LPTSolver : public Runnable, private Configuration {
 private:
  // Type declarations and template variables only
  static constexpr GInt NVARS = Particle<NDIM, P>::m_noVars;

 public:
  LPTSolver(GInt32 domainId, GInt32 noDomains) : m_domainId(domainId), m_noDomains(noDomains){};
  ~LPTSolver() override       = default;
  LPTSolver(const LPTSolver&) = delete;
  LPTSolver(LPTSolver&&)      = delete;
  auto operator=(const LPTSolver&) -> LPTSolver& = delete;
  auto operator=(LPTSolver&&) -> LPTSolver& = delete;

  void init(int argc, GChar** argv, GString config_file) override;
  void initBenchmark(int argc, GChar** argv) override;
  auto run() -> GInt override;

  /// Memory used in Kbytes
  /// \return KBytes used
  [[nodiscard]] auto mem() const -> GDouble { return m_vars.size() * 8.0 / DKBIT; }

 private:
  void init(int argc, GChar** argv);
  void initTimers();
  void allocateMemory();
  void loadConfiguration();

  void initialCondition();
  void timeStep();
  void output(const GBool forced);

  /// Configuration
  GString m_exe;
  GString m_configurationFileName;
  GBool   m_benchmark = false;
  GInt    m_capacity  = 1000;

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
  GInt                 m_timeStep;
  GDouble              m_dt;
};


#endif // LPT_SOLVER_H
