#ifndef LBM_SOLVER_H
#define LBM_SOLVER_H
#include <sfcmm_common.h>
#include "cartesiangrid.h"
#include "cell_filter.h"
#include "common/configuration.h"
#include "constants.h"
#include "interface/solver_interface.h"
#include "lbm/bnd/bnd.h"
#include "postprocess/postprocessing.h"
#include "variables.h"

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ>
class LBMSolver : public Runnable,
                  private Configuration,
                  private Postprocess<DEBUG_LEVEL, dim(LBTYPE), SolverType::LBM, LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>> {
 private:
  // Type declarations and template variables only
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
  static constexpr GInt NVAR  = noVars<LBTYPE>(EQ);

  using METH = LBMethod<LBTYPE>;

  using POST = Postprocess<DEBUG_LEVEL, dim(LBTYPE), SolverType::LBM, LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>>;

  // give postprocessing access to all data
  friend POST;
  using POST::executePostprocess;

  using VAR = LBMVariables<EQ, NDIM>;

 public:
  LBMSolver(GInt32 domainId, GInt32 noDomains) : POST(), m_domainId(domainId), m_noDomains(noDomains){};
  ~LBMSolver() override = default;

  LBMSolver(const LBMSolver&)                    = delete;
  LBMSolver(LBMSolver&&)                         = delete;
  auto operator=(const LBMSolver&) -> LBMSolver& = delete;
  auto operator=(LBMSolver&&) -> LBMSolver&      = delete;

  void init(int argc, GChar** argv, GString config_file) override;
  void initBenchmark(int argc, GChar** argv) override;
  auto run() -> GInt override;

  [[nodiscard]] auto grid() const -> const CartesianGrid<DEBUG_LEVEL, NDIM>& override {
    return *static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get());
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

  auto getCartesianGridData() -> CartesianGridData<NDIM> {
    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->getCartesianGridData();
  }


 private:
  void init(int argc, GChar** argv);
  void initTimers();
  void allocateMemory();
  void loadConfiguration();

  //  [[nodiscard]] auto grid() const -> CartesianGrid<DEBUG_LEVEL, LBMethod<LBTYPE>::m_dim>* {
  //    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get());
  //  }

  auto hasBndrySurface(const GString& surfaceName) -> GBool {
    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->hasBndrySurface(surfaceName);
  }

  auto bndrySurface(const GString& surfaceName) -> Surface<DEBUG_LEVEL, NDIM>& {
    return static_cast<CartesianGrid<DEBUG_LEVEL, NDIM>*>(m_grid.get())->bndrySurface(surfaceName);
  }


  void timeStep();
  void output(const GBool forced = false, const GString& postfix = "");
  void compareToAnalyticalResult();
  auto shiftCenter(const Point<NDIM>& center) const -> Point<NDIM>;
  void initialCondition();
  void initBndryValues();
  void writeInfo();
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
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      if(NVAR < VAR::rho()) {
        TERMM(-1, "Invalid number of variables (" + std::to_string(NVAR) + ")");
      }
      return m_vars.at(cellId * NVAR + VAR::rho());
    }
    return m_vars[cellId * NVAR + VAR::rho()];
  }

  auto inline electricPotential(const GInt cellId) -> GDouble& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      if(NVAR < VAR::electricPotential()) {
        TERMM(-1, "Invalid number of variables (" + std::to_string(NVAR) + ")");
      }
      return m_vars.at(cellId * NVAR + VAR::electricPotential());
    }
    return m_vars[cellId * NVAR + VAR::electricPotential()];
  }

  // todo: this should be mapped..(???).
  auto inline velocity(const GInt cellId) -> GDouble* {
    if(EQ == LBEquationType::Poisson) {
      TERMM(-1, "Invalid use of velocity for this equation type");
    }
    return &m_vars[cellId * NVAR + VAR::velocities()[0]];
  }

  auto inline velocity(const GInt cellId, const GInt dir) -> GDouble& {
    if(EQ == LBEquationType::Poisson) {
      TERMM(-1, "Invalid use of velocity for this equation type");
    }
    return m_vars[cellId * NVAR + VAR::velocities()[dir]];
  }

  auto inline vars(const GInt cellId, const GInt varId) -> GDouble& { return m_vars[cellId * NVAR + varId]; }

  auto inline f(const GInt cellId, const GInt dir) -> GDouble& { return m_f[cellId * NDIST + dir]; }

  auto inline feq(const GInt cellId) -> GDouble* { return &m_feq[cellId * NDIST]; }

  auto inline feq(const GInt cellId, const GInt dir) -> GDouble& { return m_feq[cellId * NDIST + dir]; }

  auto inline fold(const GInt cellId, const GInt dir) -> GDouble& {
    const GInt entryId = cellId * NDIST + dir;
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      if(entryId > static_cast<GInt>(m_fold.size())) {
        TERMM(-1, "Out of bounds!");
      }
    }
    return m_fold[entryId];
  }

  [[nodiscard]] auto inline noInternalCells() const -> GInt { return m_grid->noCells(); }
  [[nodiscard]] auto inline allCells() const -> GInt { return grid().totalSize(); }
  [[nodiscard]] auto inline maxLvl() const -> GInt { return m_grid->maxLvl(); }

  void checkDivergence();

  std::unique_ptr<GridInterface>                          m_grid;
  std::unique_ptr<LBMBndManager<DEBUG_LEVEL, LBTYPE, EQ>> m_bndManager;

  /// Configuration
  GString      m_exe;
  GString      m_configurationFileName;
  GBool        m_benchmark  = false;
  GBool        m_diverged   = false;
  LBSolverType m_solverType = LBSolverType::BGK;
  GInt         m_noSpecies  = 1;

  /// MPI
  GInt32 m_domainId  = -1;
  GInt32 m_noDomains = -1;

  /// Output
  GString                                  m_outputDir                = "out/";
  GString                                  m_solutionFileName         = "solution";
  GInt                                     m_infoInterval             = defaultInfoOutInterval;
  GInt                                     m_keepAliveMsg             = defaultKeepAlive;
  GInt                                     m_convergenceCheckInterval = defaultConvCheckInterval;
  GInt                                     m_outputSolutionInterval   = defaultSolutionInterval;
  GBool                                    m_generatePath             = true;
  std::unique_ptr<CellFilterManager<NDIM>> m_filterList               = nullptr;

  /// Variables
  std::vector<GDouble> m_f;
  std::vector<GDouble> m_feq;
  std::vector<GDouble> m_fold;
  std::vector<GDouble> m_vars;
  std::vector<GDouble> m_varsold;

  /// Reference values
  GDouble m_refLength = 1.0;
  GDouble m_refRho    = 1.0;
  GDouble m_refT      = defaultT20C;
  GDouble m_refU      = 1.0;

  GDouble m_nu              = 0;
  GDouble m_re              = 1;
  GDouble m_ma              = defaultMachNumber;
  GDouble m_latticeVelocity = 1.0;

  GDouble m_relaxTime = defaultRelaxT;
  GDouble m_omega     = 1.0 / m_relaxTime;

  GInt    m_timeStep    = 0;
  GInt    m_maxTimeStep = 0;
  GDouble m_dt          = NAN;
  GDouble m_currentTime = 0;
};

#endif // LBM_SOLVER_H
