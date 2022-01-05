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
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LPTInt], "Integration", TimeKeeper[Timers::LPTMainLoop]);
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
    RECORD_TIMER_START(TimeKeeper[Timers::LPTCalc]);
    timeStep();
    RECORD_TIMER_STOP(TimeKeeper[Timers::LPTCalc]);


    // also write an output if it's the last time step or if the solution has converged
    RECORD_TIMER_START(TimeKeeper[Timers::LPTIo]);
    output(m_timeStep == noTimesteps - 1);
    RECORD_TIMER_STOP(TimeKeeper[Timers::LPTIo]);
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
  GInt noParticles = 100;
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
  m_gravity               = {0, 10};
  m_dt                    = 0.01;
  const GInt m_maxNoSteps = 200000;

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
    // todo: make settable
    calcA<force::Model::constDensityRatioGravBuoNlinDrag, IntegrationMethod::ImplicitEuler>();
    // todo: make settable
    timeIntegration<IntegrationMethod::ImplicitEuler>();
  }

  for(GInt partId = 0; partId < m_noParticles; ++partId) {
    cerr0 << "V " << strStreamify<NDIM>(velocity(partId)).str() << std::endl;
    cerr0 << "center " << strStreamify<NDIM>(center(partId)).str() << std::endl;
  }
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
template <force::Model FM, IntegrationMethod IM>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA() {
  RECORD_TIMER_START(TimeKeeper[Timers::LPTForce]);

  GDouble       rho_a  = m_rho_a_infty;
  GDouble       nu_a   = m_nu_a_infty;
  VectorD<NDIM> velo_a = m_velo_a_infty;

  GDouble       rho_p = density(0);
  GDouble       r_p   = NAN;
  GDouble       re_p  = NAN;
  VectorD<NDIM> velo_p;
  velo_p.fill(NAN);

  for(GInt partId = 0; partId < m_noParticles; ++partId) {
    using namespace force;

    if(FM == Model::constDensityaGravBuo || FM == Model::constDensityaGravBuoStokesDrag || FM == Model::constDensityaGravBuoNlinDrag) {
      rho_p = density(partId);
    }
    if(FM == Model::constDensityRatioGravStokesDrag || FM == Model::constDensityRatioGravBuoStokesDrag
       || FM == Model::constDensityaGravBuoStokesDrag || FM == Model::constDensityRatioGravNlinDrag
       || FM == Model::constDensityRatioGravBuoNlinDrag || FM == Model::constDensityaGravBuoNlinDrag) {
      r_p = radius(partId);
    }

    if(FM == Model::constDensityRatioGravNlinDrag || FM == Model::constDensityRatioGravBuoNlinDrag
       || FM == Model::constDensityaGravBuoNlinDrag) {
      velo_p = velocity(partId);
    }

    // Switch for gravity model
    switch(FM) {
      case Model::constDensityRatioGrav:
      case Model::constDensityRatioGravStokesDrag:
      case Model::constDensityRatioGravNlinDrag:
        a(partId) = m_gravity;
        break;
      case Model::constDensityRatioGravBuo:
      case Model::constDensityRatioGravBuoStokesDrag:
      case Model::constDensityRatioGravBuoNlinDrag:
      case Model::constDensityaGravBuo:
      case Model::constDensityaGravBuoStokesDrag:
      case Model::constDensityaGravBuoNlinDrag:
        a(partId) = m_gravity * (1.0 - rho_a / rho_p);
        break;
      default:
        break;
    }
    // Switch for drag model
    // todo: reduce code repetition
    switch(FM) {
      case Model::constDensityRatioGravStokesDrag:
      case Model::constDensityRatioGravBuoStokesDrag:
      case Model::constDensityaGravBuoStokesDrag:
        DC(partId) = 18.0 * nu_a / (4.0 * r_p * r_p * rho_p);
        if(IM != IntegrationMethod::ImplicitEuler) {
          // just calculate DC since implicit methods use DC directly
          a(partId) += DC(partId) * (velo_a - velocity(partId));
        }
        break;
      case Model::constDensityRatioGravNlinDrag:
      case Model::constDensityRatioGravBuoNlinDrag:
      case Model::constDensityaGravBuoNlinDrag:
        re_p = rho_a * velo_p.norm() * 2.0 * r_p / nu_a;
        // todo: allow selecting drag model
        DC(partId) = 18.0 * nu_a / (4.0 * r_p * r_p * rho_p) * dragCoefficient<DragModel::Putnam61>(re_p);
        if(IM != IntegrationMethod::ImplicitEuler) {
          // just calculate DC since implicit methods use DC directly
          a(partId) += DC(partId) * (velo_a - velocity(partId));
        }
        break;
      default:
        DC(partId) = 0;
        break;
    }

    if(DEBUG_LEVEL >= Debug_Level::debug) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        if(isnan(a(partId)[dir]) || isinf(a(partId)[dir])) {
          TERMM(-1, "Result diverged in calcA()!");
        }
      }
    }
  }
  RECORD_TIMER_STOP(TimeKeeper[Timers::LPTForce]);
}

template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
template <IntegrationMethod IM>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::timeIntegration() {
  RECORD_TIMER_START(TimeKeeper[Timers::LPTInt]);
  VectorD<NDIM> velo_a = m_velo_a_infty;

  for(GInt partId = 0; partId < m_noParticles; ++partId) {
    const VectorD<NDIM> old_velo_p = velocity(partId);
    switch(IM) {
      case IntegrationMethod::ForwardEuler:
        velocity(partId) = old_velo_p + a(partId) * m_dt;
        center(partId)   = center(partId) + velocity(partId) * m_dt;
        break;
      case IntegrationMethod::ImplicitEuler:
        // assumes that DC and v_a is constant!
        // Derivation:
        // v1 = v0 + a(v1) * dt <->
        // v1 = v0 + g * (1.0 - rho_a / rho_p) * dt + DM * dt (gravity part is constant and replaced by G)
        // v1 = v0 + G * dt + DC * (v_a - v1) * dt <->
        // v1 = v0 + G * dt + DC * va * dt - DC * v1 * dt <->
        // v1 * (1 + DC * dt) = v0 + G * dt + DC * va * dt <->
        // v1  = (v0 + G * dt + DC * va * dt)/ (1 + DC * dt)
        velocity(partId) = (old_velo_p + a(partId) * m_dt + DC(partId) * velo_a * m_dt) / (1.0 + DC(partId) * m_dt);
        a(partId)        = (old_velo_p - velocity(partId)) / m_dt;
        center(partId)   = center(partId) + 0.5 * (velocity(partId) + old_velo_p) * m_dt;
        break;
      default:
        TERMM(-1, "Invalid Integration method!");
    }

    if(DEBUG_LEVEL >= Debug_Level::debug) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        if(isnan(velocity(partId)[dir]) || isinf(velocity(partId)[dir])) {
          cerr0 << "oldV " << strStreamify<NDIM>(old_velo_p).str() << std::endl;
          cerr0 << "a " << strStreamify<NDIM>(a(partId)).str() << std::endl;
          TERMM(-1, "Result diverged!");
        }
      }
    }
  }
  RECORD_TIMER_STOP(TimeKeeper[Timers::LPTInt]);
}

// todo: allow choice of output variables
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::output(const GBool forced, const GString& postfix) {
  if((m_timeStep > 0 && m_timeStep % m_outputSolutionInterval == 0) || forced) {
    std::vector<IOIndex>              index;
    std::vector<std::vector<GString>> values;

    std::vector<VectorD<NDIM>>             tmpCenter;
    std::array<std::vector<GDouble>, NDIM> tmpVel;
    std::vector<GDouble>                   tmpRho;
    std::vector<GDouble>                   tmpRadius;
    std::vector<GDouble>                   tmpTemperature;


    for(GInt dir = 0; dir < NDIM; ++dir) {
      tmpVel[dir].resize(m_noParticles);
    }
    tmpRho.resize(m_noParticles);
    tmpCenter.resize(m_noParticles);
    tmpRadius.resize(m_noParticles);
    tmpTemperature.resize(m_noParticles);

    for(GInt partId = 0; partId < m_noParticles; ++partId) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        tmpVel[dir][partId] = velocity(partId, dir);
      }
      tmpRho[partId]         = density(partId);
      tmpCenter[partId]      = center(partId);
      tmpRadius[partId]      = radius(partId);
      tmpTemperature[partId] = temperature(partId);
    }

    std::vector<GString> velStr = {"u", "v", "w"};
    for(GInt dir = 0; dir < NDIM; ++dir) {
      // todo: fix type
      index.emplace_back(IOIndex{velStr[dir], "float32"});
      values.emplace_back(toStringVector(tmpVel[dir], m_noParticles));
    }
    // todo: fix type
    index.emplace_back(IOIndex{"rho", "float32"});
    values.emplace_back(toStringVector(tmpRho, m_noParticles));
    // todo: fix type
    index.emplace_back(IOIndex{"radius", "float32"});
    values.emplace_back(toStringVector(tmpRadius, m_noParticles));
    // todo: fix type
    index.emplace_back(IOIndex{"temperature", "float32"});
    values.emplace_back(toStringVector(tmpTemperature, m_noParticles));

    VTK::ASCII::writePoints<NDIM>(m_outputDir + m_solutionFileName + "_" + std::to_string(m_timeStep) + postfix, m_noParticles, tmpCenter,
                                  index, values);
  }
}