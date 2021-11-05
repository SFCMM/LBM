#ifndef LBM_LBM_SOLVER_H
#define LBM_LBM_SOLVER_H
#include <sfcmm_common.h>
#include "cartesiangrid.h"
#include "interface/lbm_interface.h"
#include "interface/solver_interface.h"
#include "lbm_bnd.h"
#include "lbm_constants.h"
#include "pv.h"

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMSolver : public SolverInterface {
 public:
  LBMSolver(GInt32 domainId, GInt32 noDomains) : m_domainId(domainId), m_noDomains(noDomains){};
  ~LBMSolver() override       = default;
  LBMSolver(const LBMSolver&) = delete;
  LBMSolver(LBMSolver&&)      = delete;
  auto operator=(const LBMSolver&) -> LBMSolver& = delete;
  auto operator=(LBMSolver&&) -> LBMSolver& = delete;

  void init(int argc, GChar** argv, GString config_file) override;
  auto run() -> GInt override;

  void initBenchmark(int argc, GChar** argv) override {
    init(argc, argv);
    TERMM(-1, "Not implemented!");
  };

  [[nodiscard]] auto grid() const -> const GridInterface& override { return *m_grid; };

  template <GInt NDIM>
  [[nodiscard]] auto grid() const -> CartesianGrid<DEBUG_LEVEL, NDIM>* {
    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get());
  }

  void transferGrid(const GridInterface& grid) override;

  constexpr auto isThermal() -> GBool;

 private:
  using method = LBMethod<LBTYPE>;

  static constexpr GInt NDIM                     = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST                    = LBMethod<LBTYPE>::m_noDists;

  void init(int argc, GChar** argv);
  void initTimers();

  void setupMethod();
  void finishInit();


  void loadConfiguration();
  void timeStep();
  void initialCondition();
  void boundaryCnd();
  void propagationStep();
  void currToOldVars();
  void updateMacroscopicValues();
  void calcEquilibriumMoments();
  void collisionStep();
  [[nodiscard]] auto convergence(const GInt var) const -> GDouble;

  template <GInt NDIM>
  auto inline rho(const GInt cellId) -> GDouble& {
    if(DEBUG_LEVEL > Debug_Level::min_debug) {
      if(m_noVars < PV::rho<NDIM>()) {
        TERMM(-1, "Invalid number of variables (" + std::to_string(m_noVars) + ")");
      }
      return m_vars.at(cellId * m_noVars + PV::rho<NDIM>());
    }
    return m_vars[cellId * m_noVars + PV::rho<NDIM>()];
  }
  template <GInt NDIM>
  auto inline velocity(const GInt cellId, const GInt dir) -> GDouble& {
    return m_vars[cellId * m_noVars + PV::velocities<NDIM>()[dir]];
  }

  auto inline f(const GInt cellId, const GInt dir) -> GDouble& { return m_f[cellId * NDIST + dir]; }

  auto inline fold(const GInt cellId, const GInt dir) -> GDouble& { return m_fold[cellId * NDIST + dir]; }

  [[nodiscard]] auto inline noCells() const -> GInt { return m_grid->noCells(); }

  std::unique_ptr<GridInterface>              m_grid;
  std::unique_ptr<LBMethodInterface>          m_lbm;        // todo:implement
  std::unique_ptr<LBMBndManager<DEBUG_LEVEL>> m_bndManager; // todo:implement

  GString m_exe;
  GString m_configurationFileName;

  GInt32                m_domainId               = -1;
  GInt32                m_noDomains              = -1;

  GInt                  m_noVars                 = 1;
  GInt                  m_noSpecies              = 1;
  GInt                  m_outputInfoInterval     = defaultInfoOutInterval;
  GInt                  m_outputSolutionInterval = defaultSolutionInterval;
  GInt                  m_timeStep               = 0;

  //  LBMethodType         m_method     = LBMethodType::D2Q9;
  LBSolverType         m_solverType = LBSolverType::BGK;
//  std::vector<GDouble> m_weight;
  std::vector<GDouble> m_f;
  std::vector<GDouble> m_feq;
  std::vector<GDouble> m_fold;
  std::vector<GDouble> m_vars;
  std::vector<GDouble> m_varsold;

  // Reference values
  GDouble m_refLength = 1.0;
  GDouble m_refRho    = 1.0;
  GDouble m_refT      = defaultT20C;
  GDouble m_refU      = 1.0;

  GDouble m_nu = 0;
  GDouble m_re = 1;
  GDouble m_ma = defaultMachNumber;

  // todo: move to method impl
  GDouble m_relaxTime = defaultRelaxT;
  GDouble m_omega     = 1.0 / m_relaxTime;
};

#endif // LBM_LBM_SOLVER_H
