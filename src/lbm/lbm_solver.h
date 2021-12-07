#ifndef LBM_SOLVER_H
#define LBM_SOLVER_H
#include <sfcmm_common.h>
#include "cartesiangrid.h"
#include "common/configuration.h"
#include "interface/lbm_interface.h"
#include "interface/solver_interface.h"
#include "lbm_bnd.h"
#include "lbm_constants.h"
#include "lbm_variables.h"
#include "postprocessing.h"

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
class LBMSolver : public SolverInterface,
                  private Configuration,
                  private Postprocess<DEBUG_LEVEL, dim(LBTYPE), SolverType::LBM, LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>> {
 private:
  // Type declarations and template variables only
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
  static constexpr GInt NVARS = noVars<LBTYPE>(EQ);

  using METH = LBMethod<LBTYPE>;

  using POST = Postprocess<DEBUG_LEVEL, dim(LBTYPE), SolverType::LBM, LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>>;

  // give postprocessing access to all data
  friend POST;
  using POST::executePostprocess;

  using VAR = LBMVariables<EQ, NDIM>;

 public:
  LBMSolver(GInt32 domainId, GInt32 noDomains) : POST(), m_domainId(domainId), m_noDomains(noDomains){};
  ~LBMSolver() override       = default;
  LBMSolver(const LBMSolver&) = delete;
  LBMSolver(LBMSolver&&)      = delete;
  auto operator=(const LBMSolver&) -> LBMSolver& = delete;
  auto operator=(LBMSolver&&) -> LBMSolver& = delete;

  void init(int argc, GChar** argv, GString config_file) override;
  void initBenchmark(int argc, GChar** argv) override;
  auto run() -> GInt override;

  [[nodiscard]] auto grid() const -> const BaseCartesianGrid<DEBUG_LEVEL, NDIM>& override {
    return *static_cast<BaseCartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get());
  };
  void transferGrid(const GridInterface& grid) override;

  constexpr auto isThermal() -> GBool;

 protected:
  [[nodiscard]] auto size() const -> GInt { return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->size(); }

  [[nodiscard]] auto noLeafCells() const -> GInt { return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->noLeafCells(); }

  [[nodiscard]] auto noBndCells() const -> GInt { return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->noBndCells(); }

  [[nodiscard]] auto noChildren(const GInt cellId) const -> GInt {
    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->noChildren(cellId);
  }

  inline auto center() const -> const std::vector<Point<NDIM>>& {
    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->center();
  }

  inline auto center(const GInt cellId) const -> const Point<NDIM>& {
    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->center(cellId);
  }


  [[nodiscard]] inline auto center(const GInt cellId, const GInt dir) const -> GDouble {
    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->center(cellId, dir);
  }

  auto getCartesianGridData() -> CartesianGridData<NDIM> { return grid().getCartesianGridData(); }


 private:
  void init(int argc, GChar** argv);
  void initTimers();
  void allocateMemory();

  //  [[nodiscard]] auto grid() const -> CartesianGrid<DEBUG_LEVEL, LBMethod<LBTYPE>::m_dim>* {
  //    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get());
  //  }

  auto bndrySurface(const GString surfaceName) const -> const Surface<DEBUG_LEVEL, NDIM>& {
    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->bndrySurface(surfaceName);
  }


  void loadConfiguration();
  void timeStep();
  void output(const GBool forced = false, const GString& postfix = "");
  void compareToAnalyticalResult();
  void initialCondition();
  auto convergenceCondition() -> GBool;
  void forcing();
  void prePropBoundaryCnd();
  void boundaryCnd();
  void propagationStep();
  void currToOldVars();
  void updateMacroscopicValues();
  void calcEquilibriumMoments();
  void collisionStep();


  [[nodiscard]] auto sumAbsDiff(const GInt var) const -> GDouble;

  auto inline rho(const GInt cellId) -> GDouble& {
    if(DEBUG_LEVEL > Debug_Level::min_debug) {
      if(NVARS < VAR::rho()) {
        TERMM(-1, "Invalid number of variables (" + std::to_string(NVARS) + ")");
      }
      return m_vars.at(cellId * NVARS + VAR::rho());
    }
    return m_vars[cellId * NVARS + VAR::rho()];
  }

  auto inline electricPotential(const GInt cellId) -> GDouble& {
    if(DEBUG_LEVEL > Debug_Level::min_debug) {
      if(NVARS < VAR::electricPotential()) {
        TERMM(-1, "Invalid number of variables (" + std::to_string(NVARS) + ")");
      }
      return m_vars.at(cellId * NVARS + VAR::electricPotential());
    }
    return m_vars[cellId * NVARS + VAR::electricPotential()];
  }

  auto inline velocity(const GInt cellId, const GInt dir) -> GDouble& { return m_vars[cellId * NVARS + VAR::velocities()[dir]]; }

  auto inline vars(const GInt cellId, const GInt varId) -> GDouble& { return m_vars[cellId * NVARS + varId]; }

  auto inline f(const GInt cellId, const GInt dir) -> GDouble& { return m_f[cellId * NDIST + dir]; }

  auto inline feq(const GInt cellId, const GInt dir) -> GDouble& { return m_feq[cellId * NDIST + dir]; }

  auto inline fold(const GInt cellId, const GInt dir) -> GDouble& { return m_fold[cellId * NDIST + dir]; }

  [[nodiscard]] auto inline noCells() const -> GInt { return m_grid->noCells(); }

  std::unique_ptr<GridInterface>                      m_grid;
  std::unique_ptr<LBMethodInterface>                  m_lbm;        // todo:implement
  std::unique_ptr<LBMBndManager<DEBUG_LEVEL, LBTYPE>> m_bndManager; // todo:implement

  GString m_exe;
  GString m_configurationFileName;
  GBool   m_benchmark = false;
  GBool   m_diverged  = false;

  GInt32 m_domainId  = -1;
  GInt32 m_noDomains = -1;

  GInt m_noSpecies              = 1;
  GInt m_outputInfoInterval     = defaultInfoOutInterval;
  GInt m_outputSolutionInterval = defaultSolutionInterval;
  GInt m_timeStep               = 0;

  LBSolverType         m_solverType = LBSolverType::BGK;
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

#endif // LBM_SOLVER_H
