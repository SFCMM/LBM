#include "solver.h"
#include "globaltimers.h"

using namespace std;

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::init(int argc, GChar** argv) {
  m_exe = argv[0]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

#ifndef GRIDGEN_SINGLE_FILE_LOG
  logger.open("lpt_log" + std::to_string(m_domainId), false, argc, argv, MPI_COMM_WORLD);
#else
  if(DEBUG_LEVEL < Debug_Level::max_debug) {
    logger.open("lpt_log", true, argc, argv, MPI_COMM_WORLD);
  } else {
    logger.open("lpt_log", false, argc, argv, MPI_COMM_WORLD);
  }
#endif
  logger.setMinFlushSize(LOG_MIN_FLUSH_SIZE);

  initTimers();
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::init(int argc, GChar** argv, GString config_file) {
  setConfiguration(config_file);
  init(argc, argv);
  logger << NDIM << "D LPT Solver started ||>" << endl;
  cout << NDIM << "D LPT Solver started ||>" << endl;
  Configuration::load("solver");
  //  POST::setConfAccessor(Configuration::getAccessor("postprocessing"));
  RECORD_TIMER_STOP(TimeKeeper[Timers::LPTInit]);
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::initBenchmark(int argc, GChar** argv) {
  m_benchmark = true;
  init(argc, argv);

  logger << "Setting up benchmarking!" << endl;

  TERMM(-1, "Not implemented!");
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::initTimers() {
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LPTSolverTotal], "Total run time of the LPT Solver.", TimeKeeper[Timers::timertotal]);
  RECORD_TIMER_START(TimeKeeper[Timers::LPTSolverTotal]);

  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LPTInit], "Initialization of the LPT solver", TimeKeeper[Timers::LPTSolverTotal]);
  RECORD_TIMER_START(TimeKeeper[Timers::LPTInit]);

  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LPTMainLoop], "Main Loop of the LPT solver!", TimeKeeper[Timers::LPTSolverTotal]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LPTCalc], "Computation", TimeKeeper[Timers::LPTMainLoop]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LPTColl], "Collision", TimeKeeper[Timers::LPTCalc]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LPTForce], "Computing forces", TimeKeeper[Timers::LPTCalc]);


  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LPTPost], "Postprocessing", TimeKeeper[Timers::LPTMainLoop]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LPTIo], "IO", TimeKeeper[Timers::LPTMainLoop]);
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::loadConfiguration() {
  m_outputDir              = opt_config_value<GString>("output_dir", m_outputDir);
  m_outputInfoInterval     = opt_config_value<GInt>("info_interval", m_outputInfoInterval);
  m_outputSolutionInterval = opt_config_value<GInt>("solution_interval", m_outputSolutionInterval);
  m_solutionFileName       = opt_config_value<GString>("solution_filename", m_solutionFileName);
  m_initialCondition       = getLPTInitCond(opt_config_value<GString>("initialcondition", "none"));
  // todo: set by configuration
  m_capacity = 1000;

  if(!isPath(m_outputDir, m_generatePath)) {
    TERMM(-1, "Invalid output directory set! (value: " + m_outputDir + ")");
  }
  // fix path
  if(m_outputDir.back() != '/') {
    m_outputDir += '/';
  }
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::allocateMemory() {
  // allocate memory
  m_vars.resize(m_capacity * NVARS);
  // fill variables with NAN to crash when using unset variables
  fill(m_vars.begin(), m_vars.end(), NAN);
  cerr0 << "allocated memory: " << mem() << "KB" << std::endl;
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
auto LPTSolver<DEBUG_LEVEL, NDIM, P>::run() -> GInt {
  RECORD_TIMER_START(TimeKeeper[Timers::LPTInit]);
  loadConfiguration();
  allocateMemory();
  //  POST::init();
  initialCondition();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LPTInit]);

  RECORD_TIMER_START(TimeKeeper[Timers::LPTMainLoop]);
  // todo: make settable
  const GInt noTimesteps = required_config_value<GInt>("maxSteps");

  //  executePostprocess(pp::HOOK::ATSTART);

  for(m_timeStep = 0; m_timeStep < noTimesteps; ++m_timeStep) {
    timeStep();

    // also write an output if it's the last time step or if the solution has converged
    output(m_timeStep == noTimesteps - 1);
  }

  //  if(has_config_value("analyticalSolution")) {
  //    compareToAnalyticalResult();
  //  }

  //  executePostprocess(pp::HOOK::ATEND);

  logger << "LPT Solver finished <||" << endl;
  cout << "LPT Solver finished <||" << endl;
  unusedConfigValues();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LPTMainLoop]);
  return 0;
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::initialCondition() {
  switch(m_initialCondition) {
    case LPTInitCond::randomvol_pos:
      init_randomVolPos();
      break;
    case LPTInitCond::none:
      return;
    default:
      TERMM(-1, "Invalid initial condition");
  }
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::init_randomVolPos() {
  // todo: make settable
  GInt noParticles = 10;
  // todo: make settable
  GString m_volType = "box";
  // todo: make settable
  VectorD<NDIM> cornerA = {0, 0};
  // todo: make settable
  VectorD<NDIM> cornerB = {10, 1};

  cerr0 << "Placing " << noParticles << " particles inside " << m_volType << " at random position" << std::endl;

  for(GInt partId = 0; partId < noParticles; ++partId) {
    center(partId) = randomPos<NDIM>(cornerA, cornerB);
    velocity(partId).fill(0);
  }
}


template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::timeStep() {
  // todo: load from config
  VectorD<NDIM> m_gravity     = {0, 10};
  GInt          m_noParticles = 10;
  GDouble       m_dt          = 0.1;
  GInt          m_maxNoSteps  = 100;

  for(m_timeStep = 0; m_timeStep < m_maxNoSteps; ++m_timeStep) {
    for(GInt partId = 0; partId < m_noParticles; ++partId) {
      velocity(partId) = m_gravity * m_dt;
      center(partId)   = center(partId) + velocity(partId) * m_dt;
    }
  }

  for(GInt partId = 0; partId < m_noParticles; ++partId) {
    cerr0 << "V " << strStreamify<NDIM>(velocity(partId)).str() << std::endl;
    cerr0 << "center " << strStreamify<NDIM>(center(partId)).str() << std::endl;
  }
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::output(const GBool forced) {}