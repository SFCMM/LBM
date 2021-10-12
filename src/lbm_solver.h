#ifndef LBM_LBM_SOLVER_H
#define LBM_LBM_SOLVER_H
#include <sfcmm_common.h>
#include "cartesiangrid.h"
#include "interface/app_interface.h"
#include "lbm_constants.h"

template <Debug_Level DEBUG_LEVEL>
class LBMSolver : public AppInterface {
 public:
  LBMSolver(GInt32 domainId, GInt32 noDomains) : m_domainId(domainId), m_noDomains(noDomains){};
  ~LBMSolver() override = default;

  void init(int argc, GChar** argv, GString config_file) override;


  void initBenchmark(int argc, GChar** argv) override {
    init(argc, argv);
    TERMM(-1, "Not implemented!");
  };

  auto               run() -> GInt override;
  [[nodiscard]] auto grid() const -> const GridInterface& override { return *m_grid; };

  void transferGrid(const GridInterface& grid) override;

  constexpr auto isThermal() -> GBool;

 private:
  void init(int argc, GChar** argv);
  void initTimers();

  std::unique_ptr<GridInterface> m_grid;

  GString m_exe;
  GString m_configurationFileName;

  GInt32 m_domainId  = -1;
  GInt32 m_noDomains = -1;
  GInt   m_dim       = 0;
  GInt   m_noVars    = 1;
  GInt   m_noSpecies = 1;

  LBMethodType m_method     = LBMethodType::D2Q9;
  LBSolverType m_solverType = LBSolverType::BGK;

  std::vector<GDouble> m_f;
  std::vector<GDouble> m_feq;
  std::vector<GDouble> m_fold;
  std::vector<GDouble> m_vars;
  std::vector<GDouble> m_varsold;

  void loadConfiguration();
  void setupMethod();
  void finishInit();

  template <GInt NDIM>
  void timeStep();
  void initialCondition();
};

#endif // LBM_LBM_SOLVER_H
