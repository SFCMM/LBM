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
auto LBMSolver<DEBUG_LEVEL>::run() -> GInt {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMMainLoop]);

  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMMainLoop]);
  return 0;
}

template class LBMSolver<Debug_Level::no_debug>;
template class LBMSolver<Debug_Level::min_debug>;
template class LBMSolver<Debug_Level::debug>;
template class LBMSolver<Debug_Level::more_debug>;
template class LBMSolver<Debug_Level::max_debug>;