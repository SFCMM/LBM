#ifndef LBM_SOLVEREXE_H
#define LBM_SOLVEREXE_H

#include <sfcmm_common.h>
#include "lbm_solver.h"

// encapsulation to specialize to LBMethodeType
template <Debug_Level DEBUG_LEVEL>
class LBMSolverExecutor : public SolverInterface {
 public:
  LBMSolverExecutor(const GInt32 domainId, const GInt32 noDomains) : m_domainId(domainId), m_noDomains(noDomains){};

  void init(int argc, GChar** argv, GString config_file) override {
    static constexpr GInt nodims = 2;
    static constexpr GInt ndist  = 9;

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
    LBEquation equationType = getLBEquationType(equation);

    switch(modelType) {
      case LBMethodType::D1Q3:
        switch(equationType) {
          case LBEquation::Navier_Stokes:
            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D1Q3, LBEquation::Navier_Stokes>>(m_domainId, m_noDomains);
            break;
          case LBEquation::Poisson:
            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D1Q3, LBEquation::Poisson>>(m_domainId, m_noDomains);
            break;
          case LBEquation::Navier_Stokes_Poisson:
            m_lbmSolver =
                std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D1Q3, LBEquation::Navier_Stokes_Poisson>>(m_domainId, m_noDomains);
            break;
          default:
            TERMM(-1, "Unsupported equation type");
        }
        break;
      case LBMethodType::D2Q5:
        switch(equationType) {
          case LBEquation::Navier_Stokes:
            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D2Q5, LBEquation::Navier_Stokes>>(m_domainId, m_noDomains);
            break;
          case LBEquation::Poisson:
            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D2Q5, LBEquation::Poisson>>(m_domainId, m_noDomains);
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
          case LBEquation::Navier_Stokes:
            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D2Q9, LBEquation::Navier_Stokes>>(m_domainId, m_noDomains);
            break;
          case LBEquation::Poisson:
            m_lbmSolver = std::make_unique<LBMSolver<DEBUG_LEVEL, LBMethodType::D2Q9, LBEquation::Poisson>>(m_domainId, m_noDomains);
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
  void initBenchmark(int argc, GChar** argv) override { m_lbmSolver->initBenchmark(argc, argv); }

  auto run() -> GInt override { return m_lbmSolver->run(); }

  [[nodiscard]] auto grid() const -> const GridInterface& override { return m_lbmSolver->grid(); }
  void               transferGrid(const GridInterface& grid) override { m_lbmSolver->transferGrid(grid); }

 private:
  GInt32 m_domainId  = -1;
  GInt32 m_noDomains = -1;

  std::unique_ptr<SolverInterface> m_lbmSolver;
};
#endif // LBM_SOLVEREXE_H
