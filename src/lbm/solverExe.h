#ifndef LBM_SOLVEREXE_H
#define LBM_SOLVEREXE_H

#include <sfcmm_common.h>
#include "solver.h"

// encapsulation to specialize to LBMethodeType
template <Debug_Level DEBUG_LEVEL>
class LBMSolverExecutor : public Runnable {
 public:
  LBMSolverExecutor(const GInt32 domainId, const GInt32 noDomains) : m_domainId(domainId), m_noDomains(noDomains){};

  /// Init the LBM solver and create the correct templated version of the solver.
  /// \param argc CMD
  /// \param argv CMD
  /// \param config_file configuration file to be used
  void init(int argc, GChar** argv, GString config_file) override {
    // default model
    GString model = "D2Q9";

    // default equation
    GString equation = "navierstokes";

    if(MPI::isRoot()) {
      std::ifstream configFileStream(config_file);
      json          config;
      configFileStream >> config;
      // determine which model to use and equation to solve
      model    = config::opt_config_value(config["solver"], "model", model);
      equation = config::opt_config_value(config["solver"], "equation", equation);
    }
    LBMethodType modelType = getLBMethodType(model);
    // todo: communicate the model
    LBEquationType equationType = getLBEquationType(equation);

    switch(modelType) {
      case LBMethodType::D1Q3:
        switch(equationType) {
            //          case LBEquation::Navier_Stokes:
            //            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D1Q3, LBEquation::Navier_Stokes>>(m_domainId,
            //            m_noDomains); break;
          case LBEquationType::Poisson:
            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D1Q3, LBEquationType::Poisson>>(m_domainId, m_noDomains);
            break;
            //          case LBEquation::Navier_Stokes_Poisson:
            //            m_lbmSolver =
            //                std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D1Q3, LBEquation::Navier_Stokes_Poisson>>(m_domainId,
            //                m_noDomains);
            //            break;
          default:
            TERMM(-1, "Unsupported equation type");
        }
        break;
      case LBMethodType::D2Q5:
        switch(equationType) {
            //          case LBEquation::Navier_Stokes:
            //            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D2Q5, LBEquation::Navier_Stokes>>(m_domainId,
            //            m_noDomains); break;
          case LBEquationType::Poisson:
            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D2Q5, LBEquationType::Poisson>>(m_domainId, m_noDomains);
            break;
            //          case LBEquation::Navier_Stokes_Poisson:
            //            m_lbmSolver =
            //                std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D2Q5, LBEquation::Navier_Stokes_Poisson>>(m_domainId,
            //                m_noDomains);
            //            break;
          default:
            TERMM(-1, "Unsupported equation type");
        }
        break;
      case LBMethodType::D2Q9:
        switch(equationType) {
          case LBEquationType::Navier_Stokes:
            m_lbmSolver =
                std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D2Q9, LBEquationType::Navier_Stokes>>(m_domainId, m_noDomains);
            break;
          case LBEquationType::Poisson:
            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D2Q9, LBEquationType::Poisson>>(m_domainId, m_noDomains);
            break;
            //          case LBEquation::Navier_Stokes_Poisson:
            //            m_lbmSolver =
            //                std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D2Q9, LBEquation::Navier_Stokes_Poisson>>(m_domainId,
            //                m_noDomains);
            //            break;
          default:
            TERMM(-1, "Unsupported equation type");
        }
        break;
      default:
        TERMM(-1, "Unsupported model");
    }
    m_lbmSolver->init(argc, argv, config_file);
  }

  /// Init a benchmark run for the LBM solver
  /// \param argc CMD
  /// \param argv CMD
  void initBenchmark(int argc, GChar** argv) override { m_lbmSolver->initBenchmark(argc, argv); }

  /// Run the LBM solver
  /// \return Errorcode
  auto run() -> GInt override { return m_lbmSolver->run(); }

  /// Access to the grid
  /// \return Return grid access
  [[nodiscard]] auto grid() const -> const GridInterface& override { return m_lbmSolver->grid(); }

  /// Transfer a grid which is already in memory
  /// \param grid Grid to be transferred
  void transferGrid(const GridInterface& grid) override { m_lbmSolver->transferGrid(grid); }

 private:
  GInt32 m_domainId  = -1;
  GInt32 m_noDomains = -1;

  std::unique_ptr<Runnable> m_lbmSolver;
};
#endif // LBM_SOLVEREXE_H
