#include "solver.h"
#include "common/sphere.h"
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
  // todo: make configable
  const GDouble initVelo = 0;
  // todo: make configable
  const GDouble initDensity = 2;
  // todo: make configable
  const GDouble initRadius = 0.01;

  cerr0 << "Placing " << noParticles << " particles inside " << m_volType << " at random position" << std::endl;

  for(GInt partId = 0; partId < noParticles; ++partId) {
    center(partId) = randomPos<NDIM>(cornerA, cornerB);
    velocity(partId).fill(initVelo);
    density(partId) = initDensity;
    radius(partId)  = initRadius;
    volume(partId)  = sphere::volumeR(initRadius);
    m_noParticles++;
  }
}


template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::timeStep() {
  // todo: load from config
  m_gravity                  = {0, 10};
  const GDouble m_dt         = 0.1;
  const GInt    m_maxNoSteps = 200000;

  // todo: load from config
  const GDouble rho_a        = 1;
  m_rho_a_infty              = 1;
  const GDouble nu_a         = 1E-5;
  m_nu_a_infty               = 1E-5;
  const VectorD<NDIM> velo_a = {0, 0};
  m_velo_a_infty             = {0, 0};

  for(GInt partId = 0; partId < m_noParticles; ++partId) {
    cerr0 << "V " << strStreamify<NDIM>(velocity(partId)).str() << std::endl;
    cerr0 << "center " << strStreamify<NDIM>(center(partId)).str() << std::endl;
  }

  for(m_timeStep = 0; m_timeStep < m_maxNoSteps; ++m_timeStep) {
    // 1. Step update acceleration
    calcA<force::Model::constDensityRatioGravBuoStokesDrag>();
    for(GInt partId = 0; partId < m_noParticles; ++partId) {
      const VectorD<NDIM> old_velo_p = velocity(partId);
      //      const VectorD<NDIM> old_rel_velo = velo_a - old_velo_p;
      //      const GDouble reynoldsNumber = rho_a * old_rel_velo.norm() * 2.0 * r_p / nu_a;
      //      const GDouble DC             = 18.0 * nu_a / (4.0 * r_p * r_p * rho_p);
      //      const GDouble mass = volume(partId) * rho_p;

      // Note: all forces are based on a=forces/mass
      //      a(partId) = 18.0 * nu_a / (4.0 * r_p * r_p * rho_p) * old_rel_velo // stokes drag force
      //                  + m_gravity * (1.0 - rho_a / rho_p);                   // gravity + buoyancy
      velocity(partId) = old_velo_p + a(partId) * m_dt;
      center(partId)   = center(partId) + velocity(partId) * m_dt;
      if(DEBUG_LEVEL >= Debug_Level::debug) {
        for(GInt dir = 0; dir < NDIM; ++dir) {
          if(isnan(velocity(partId)[dir]) || isinf(velocity(partId)[dir])) {
            TERMM(-1, "Result diverged!");
          }
        }
      }
    }
  }

  for(GInt partId = 0; partId < m_noParticles; ++partId) {
    cerr0 << "V " << strStreamify<NDIM>(velocity(partId)).str() << std::endl;
    cerr0 << "center " << strStreamify<NDIM>(center(partId)).str() << std::endl;
  }
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
template <force::Model FM>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA() {
  GDouble       rho_a  = m_rho_a_infty;
  GDouble       nu_a   = m_nu_a_infty;
  VectorD<NDIM> velo_a = m_velo_a_infty;

  GDouble rho_p = density(0);

  for(GInt partId = 0; partId < m_noParticles; ++partId) {
    using namespace force;

    if(FM == Model::constDensityaGravBuo || FM == Model::constDensityaGravBuoStokesDrag) {
      rho_p = density(partId);
    }

    // Switch for gravity model
    switch(FM) {
      case Model::constDensityRatioGrav:
      case Model::constDensityRatioGravStokesDrag:
        a(partId) = m_gravity;
        break;
      case Model::constDensityRatioGravBuo:
      case Model::constDensityRatioGravBuoStokesDrag:
      case Model::constDensityaGravBuo:
      case Model::constDensityaGravBuoStokesDrag:
        a(partId) = m_gravity * (1.0 - rho_a / rho_p);
        break;
      default:
        break;
    }
    // Switch for drag model
    switch(FM) {
      case Model::constDensityRatioGravStokesDrag:
      case Model::constDensityRatioGravBuoStokesDrag:
      case Model::constDensityaGravBuoStokesDrag:
        a(partId) += 18.0 * nu_a / (4.0 * radius(partId) * radius(partId) * rho_p) * (velo_a - velocity(partId));
        break;
      default:
        break;
    }
  }
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::output(const GBool forced) {}