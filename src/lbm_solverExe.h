#ifndef LBM_LBM_SOLVEREXE_H
#define LBM_LBM_SOLVEREXE_H

#include <sfcmm_common.h>
#include "lbm_solver.h"

// encapsulation to specialize to LBMethodeType
template <Debug_Level DEBUG_LEVEL>
class LBMSolverExecutor : public SolverInterface {
 public:
  LBMSolverExecutor(const GInt32 domainId, const GInt32 noDomains)
    : m_domainId(domainId),
      m_noDomains(noDomains){
    // todo: make settable
    static constexpr GInt nodims = 2;
    static constexpr GInt ndist = 9;

    m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, getLBMethodType(nodims, ndist)>>(domainId, noDomains);
  };

  void init(int argc, GChar** argv, GString config_file) override { m_lbmSolver->init(argc, argv, config_file); }
  void initBenchmark(int argc, GChar** argv) override { m_lbmSolver->initBenchmark(argc, argv); }

  auto run() -> GInt override { return m_lbmSolver->run(); }

  [[nodiscard]] auto grid() const -> const GridInterface& override { return m_lbmSolver->grid(); }
  void               transferGrid(const GridInterface& grid) override { m_lbmSolver->transferGrid(grid); }

 private:
  GInt32 m_domainId  = -1;
  GInt32 m_noDomains = -1;

  std::unique_ptr<SolverInterface> m_lbmSolver;
};
#endif // LBM_LBM_SOLVEREXE_H
