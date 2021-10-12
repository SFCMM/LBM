#include "lbm_solver.h"

using namespace std;

template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::init(int argc, GChar** argv) {
  m_exe = argv[0]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

#ifndef GRIDGEN_SINGLE_FILE_LOG
  logger.open("lbm_log" + std::to_string(m_domainId), false, argc, argv, MPI_COMM_WORLD);
#else
  if(DEBUG_LEVEL < Debug_Level::more_debug) {
    logger.open("lbm_log", true, argc, argv, MPI_COMM_WORLD);
  } else {
    logger.open("lbm_log", false, argc, argv, MPI_COMM_WORLD);
  }
#endif
  logger.setMinFlushSize(LOG_MIN_FLUSH_SIZE);

  initTimers();
}

template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::init(int argc, GChar** argv, GString config_file) {
  m_configurationFileName = std::move(config_file);
  init(argc, argv);
  logger << "LBM Solver started ||>" << endl;
  cout << "LBM Solver started ||>" << endl;
  loadConfiguration();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMInit]);
}

template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::initTimers() {
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMSolverTotal], "Total run time of the LBM Solver.", TimeKeeper[Timers::timertotal]);
  RECORD_TIMER_START(TimeKeeper[Timers::LBMSolverTotal]);

  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMInit], "Initialization of the LBM solver!", TimeKeeper[Timers::LBMSolverTotal]);
  RECORD_TIMER_START(TimeKeeper[Timers::LBMInit]);

  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMMainLoop], "Main Loop of the LBM solver!", TimeKeeper[Timers::LBMSolverTotal]);
}



template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::loadConfiguration() {
  m_solverType = LBSolverType::BGK;  // todo: load from file
  m_method     = LBMethodType::D2Q9; // todo: load from file
  m_noVars = m_dim + 1 + static_cast<GInt>(isThermal());
}



template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::finishInit() {
  //allocate memory
  m_f.resize(grid().size() * LBMethod::noDists.at(static_cast<GInt>(m_method)));
  m_feq.resize(grid().size() * LBMethod::noDists.at(static_cast<GInt>(m_method)));
  m_fold.resize(grid().size() * LBMethod::noDists.at(static_cast<GInt>(m_method)));
  m_vars.resize(grid().size() * m_noVars);
  m_varsold.resize(grid().size() * m_noVars);
  setupMethod();
}


template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::setupMethod() {}


template <Debug_Level DEBUG_LEVEL>
auto LBMSolver<DEBUG_LEVEL>::run() -> GInt {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMMainLoop]);
  initialCondition();
  const GInt noTimesteps = 100;
  for(GInt ts = 0; ts < noTimesteps; ++ts) {
    switch(m_dim) {
      case 1:
        timeStep<1>();
        break;
      case 2:
        timeStep<2>();
        break;
      case 3:
        timeStep<3>();
        break;
      case 4:
        timeStep<4>();
        break;
      default:
        TERMM(-1, "Invalid number of dimensions 1-4.");
    }
  }

  logger << "LBM Solver finished <||" << endl;
  cout << "LBM Solver finished <||" << endl;
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMMainLoop]);
  return 0;
}

template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::initialCondition() {}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::timeStep() {

}

template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::transferGrid(const GridInterface& grid) {
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

  finishInit();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMInit]);
}

template <Debug_Level DEBUG_LEVEL>
constexpr auto LBMSolver<DEBUG_LEVEL>::isThermal() -> GBool{
  return false;
}



template class LBMSolver<Debug_Level::no_debug>;
template class LBMSolver<Debug_Level::min_debug>;
template class LBMSolver<Debug_Level::debug>;
template class LBMSolver<Debug_Level::more_debug>;
template class LBMSolver<Debug_Level::max_debug>;