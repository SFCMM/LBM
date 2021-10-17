#ifndef LBM_LBM_SOLVER_H
#define LBM_LBM_SOLVER_H
#include <sfcmm_common.h>
#include "cartesiangrid.h"
#include "interface/app_interface.h"
#include "interface/lbm_interface.h"
#include "lbm_constants.h"
#include "pv.h"

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

  void loadConfiguration();
  void setupMethod();
  void finishInit();

  template <GInt NDIM>
  void timeStep();
  template <GInt NDIM>
  void initialCondition();
  template <GInt NDIM>
  void boundaryCnd();
  template <GInt NDIM>
  void propagationStep();
  template <GInt NDIM>
  void currToOldVars();
  template <GInt NDIM>
  void updateMacroscopicValues();
  template <GInt NDIM>
  void calcEquilibriumMoments();
  template <GInt NDIM>
  void collisionStep();

  template <GInt NDIM>
  auto inline rho(const GInt cellId) -> GDouble& {
    return m_vars[cellId * m_noVars + PV::rho<NDIM>()];
  }
  template <GInt NDIM>
  auto inline velocity(const GInt cellId, const GInt dir) -> GDouble& {
    return m_vars[cellId * m_noVars + PV::velocities<NDIM>()[dir]];
  }


  [[nodiscard]] auto inline noCells() const -> GInt { return m_grid->noCells(); }

  std::unique_ptr<GridInterface>     m_grid;
  std::unique_ptr<LBMethodInterface> m_lbm;

  GString m_exe;
  GString m_configurationFileName;

  GInt32 m_domainId  = -1;
  GInt32 m_noDomains = -1;
  GInt   m_dim       = 0;
  GInt   m_noVars    = 1;
  GInt   m_noSpecies = 1;
  GInt   m_noDists = 0;

  LBMethodType m_method     = LBMethodType::D2Q9;
  LBSolverType m_solverType = LBSolverType::BGK;
  std::vector<GDouble> m_weight;

  std::vector<GDouble> m_f;
  std::vector<GDouble> m_feq;
  std::vector<GDouble> m_fold;
  std::vector<GDouble> m_vars;
  std::vector<GDouble> m_varsold;

  // Reference values
  GDouble m_refLength = 1.0;
  GDouble m_refRho = 1.0;
  GDouble m_refT = 293.15;
  GDouble m_refU = 1.0;

  GDouble m_nu = 0;
  GDouble m_re = 1;
  GDouble m_ma = 0.2;

  //todo: move to method impl
  GDouble m_relaxTime = 0.9;
  GDouble m_omega = 1.0/m_relaxTime;
};

#endif // LBM_LBM_SOLVER_H
