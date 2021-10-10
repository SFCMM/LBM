#ifndef LBM_LBM_SOLVER_H
#define LBM_LBM_SOLVER_H
#include <sfcmm_common.h>
#include "cartesiangrid.h"
#include "interface/app_interface.h"
//#include "interface/grid_interface.h"

template <Debug_Level DEBUG_LEVEL>
class LBMSolver : public AppInterface {
 public:
  LBMSolver(GInt32 domainId, GInt32 noDomains) : m_domainId(domainId), m_noDomains(noDomains) {

  };
  ~LBMSolver() override = default;

  void init(int argc, GChar** argv);
  void init(int argc, GChar** argv, GString config_file) override;


  void initBenchmark(int argc, GChar** argv) override {
    init(argc, argv);
    TERMM(-1, "Not implemented!");
  };

  auto               run() -> GInt override;
  [[nodiscard]] auto grid() const -> const GridInterface& override { return *m_grid; };

  void transferGrid(const GridInterface& grid) override {
    RECORD_TIMER_START(TimeKeeper[Timers::LBMInit]);
    cerr0 << "Transferring Grid to LBM solver" << std::endl;
    logger << "Transferring Grid to LBM solver" << std::endl;

    m_dim = grid.dim();

    switch(m_dim) {
      case 1:
        m_grid = std::make_unique<CartesianGrid<DEBUG_LEVEL, 1>>();
        break;
      case 2:
        m_grid = std::make_unique<CartesianGrid<DEBUG_LEVEL, 2>>();
        break;
      case 3:
        m_grid = std::make_unique<CartesianGrid<DEBUG_LEVEL, 3>>();
        break;
      default:
        TERMM(-1, "Only dimensions 1,2 and 3 are supported.");
    }

    switch(m_dim) {
      case 1:
        static_cast<CartesianGrid<DEBUG_LEVEL, 1>*>(m_grid.get())
            ->loadGridInplace(*static_cast<const CartesianGridGen<DEBUG_LEVEL, 1>*>(static_cast<const void*>(&grid)));
        break;
      case 2:
        static_cast<CartesianGrid<DEBUG_LEVEL, 2>*>(m_grid.get())
            ->loadGridInplace(*static_cast<const CartesianGridGen<DEBUG_LEVEL, 2>*>(static_cast<const void*>(&grid)));
        break;
      case 3:
        static_cast<CartesianGrid<DEBUG_LEVEL, 3>*>(m_grid.get())
            ->loadGridInplace(*static_cast<const CartesianGridGen<DEBUG_LEVEL, 3>*>(static_cast<const void*>(&grid)));
        break;
      default:
        TERMM(-1, "Only dimensions 1,2 and 3 are supported.");
    }
    RECORD_TIMER_STOP(TimeKeeper[Timers::LBMInit]);
  };

 private:
  void initTimers();

  std::unique_ptr<GridInterface> m_grid;

  GString m_exe;
  GString m_configurationFileName;

  GInt32 m_domainId  = -1;
  GInt32 m_noDomains = -1;
  GInt   m_dim       = 0;

};

#endif // LBM_LBM_SOLVER_H
