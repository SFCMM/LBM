#ifndef LPT_SOLVEREXE_H
#define LPT_SOLVEREXE_H

#include <sfcmm_common.h>
#include "solver.h"

// encapsulation to specialize LPT
template <Debug_Level DEBUG_LEVEL>
class LPTSolverExecutor : public Runnable {
 public:
  LPTSolverExecutor(const GInt32 domainId, const GInt32 noDomains) : m_domainId(domainId), m_noDomains(noDomains){};

  /// Init the LBM solver and create the correct templated version of the solver.
  /// \param argc CMD
  /// \param argv CMD
  /// \param config_file configuration file to be used
  void init(int argc, GChar** argv, GString config_file) override {
    //    if(MPI::isRoot()) {
    //      std::ifstream configFileStream(config_file);
    //      json          config;
    //      configFileStream >> config;
    //      // determine which model to use and equation to solve
    //      model    = config::opt_config_value(config["solver"], "model", model);
    //      equation = config::opt_config_value(config["solver"], "equation", equation);
    //    }
    //    // todo: communicate the model

    m_lptSolver = std::make_unique<LPTSolver<DEBUG_LEVEL, 2, LPTType::Normal>>(m_domainId, m_noDomains);

    m_lptSolver->init(argc, argv, config_file);
  }

  /// Init a benchmark run for the LPT solver
  /// \param argc CMD
  /// \param argv CMD
  void initBenchmark(int argc, GChar** argv) override { m_lptSolver->initBenchmark(argc, argv); }

  /// Run the LPT solver
  /// \return Errorcode
  auto run() -> GInt override { return m_lptSolver->run(); }

  /// Access to the grid
  /// \return Return grid access
  [[nodiscard]] auto grid() const -> const GridInterface& override { return m_lptSolver->grid(); }

  /// Transfer a grid which is already in memory
  /// \param grid Grid to be transferred
  void transferGrid(const GridInterface& grid) override { m_lptSolver->transferGrid(grid); }

 private:
  GInt32 m_domainId  = -1;
  GInt32 m_noDomains = -1;

  std::unique_ptr<Runnable> m_lptSolver;
};


#endif // LPT_SOLVEREXE_H
