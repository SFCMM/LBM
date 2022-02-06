#include "solver.h"
#include "analytical_solutions.h"
#include "common/sphere.h"
#include "globaltimers.h"

using namespace std;

/// Common init for both benchmark and normal run.
/// \tparam DEBUG_LEVEL Debug Mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
/// \param argc cmd argument count
/// \param argv cmd argument list
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

/// Entry point for normal simulations
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
/// \param argc cmd arguments count
/// \param argv cmd argument list
/// \param config_file config file name
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

/// Entry point for benchmark runs
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
/// \param argc cmd arguments count
/// \param argv cmd argument list
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::initBenchmark(int argc, GChar** argv) {
  m_benchmark = true;
  init(argc, argv);

  logger << "Setting up benchmarking!" << endl;

  TERMM(-1, "Not implemented!");
}

/// Init the timers for this solver
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
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

/// Load the configuration from the configuration file
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::loadConfiguration() {
  m_outputDir              = opt_config_value<GString>("output_dir", m_outputDir);
  m_outputInfoInterval     = opt_config_value<GInt>("info_interval", m_outputInfoInterval);
  m_outputSolutionInterval = opt_config_value<GInt>("solution_interval", m_outputSolutionInterval);
  m_solutionFileName       = opt_config_value<GString>("solution_filename", m_solutionFileName);
  m_initialCondition       = initCond(opt_config_value<GString>("initialcondition", "none"));

  m_capacity = opt_config_value<GInt>("capacity", m_capacity);

  m_gravity     = required_config_value<NDIM>("gravity");
  m_rho_a_infty = opt_config_value<GDouble>("rho_infty", m_rho_a_infty);
  m_nu_a_infty  = opt_config_value<GDouble>("nu_infty", m_nu_a_infty);

  std::vector<GDouble> tmp;
  tmp.resize(NDIM);
  fill_n(tmp.begin(), NDIM, 0);
  m_velo_a_infty = eigenutil::unpack<NDIM>(opt_config_value<std::vector<GDouble>>("v_infty", tmp));

  m_dt         = required_config_value<GDouble>("dt");
  m_maxNoSteps = required_config_value<GInt>("maxSteps");

  GString generationMethodName = "none";
  m_generationMethod           = generationMethod(opt_config_value("generation_method", generationMethodName));

  // if some is set init the method
  if(m_generationMethod != GenerationMethod::None) {
    initGenerationMethod();
  }

  if(!isPath(m_outputDir, m_generatePath)) {
    TERMM(-1, "Invalid output directory set! (value: " + m_outputDir + ")");
  }
  // fix path
  if(m_outputDir.back() != '/') {
    m_outputDir += '/';
  }

  setMethods();

  cerr0 << "<<<<<<<<<<<<>>>>>>>>>>>>>" << std::endl;
  cerr0 << "Generation method " << GenerationMethodName[static_cast<GInt>(m_generationMethod)] << std::endl;
  cerr0 << "Output interval " << m_outputSolutionInterval << std::endl;
  cerr0 << "+++++++++++++++++++++++++" << std::endl;
}

/// Init the generation methods by loading there configuraition
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::initGenerationMethod() {
  if(m_generationMethod != GenerationMethod::None) {
    switch(m_generationMethod) {
      case GenerationMethod::ConstantRate:
        break;
      case GenerationMethod::InjectionModel:
        m_injectors.emplace_back(Injector<NDIM>(getObject("injection_model")));
        break;
      default:
        TERMM(-1, "Invalid generation method selected.");
    }
  }
}

/// Set the time integration method
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::setMethods() {
  GString integrationName = opt_config_value("integrationMethod", static_cast<GString>("ForwardEuler"));

  switch(integrationMethod(integrationName)) {
    case IntegrationMethod::ForwardEuler:
      m_timeIntegration = &LPTSolver<DEBUG_LEVEL, NDIM, P>::timeIntegration<IntegrationMethod::ForwardEuler>;
      setForceModel<IntegrationMethod::ForwardEuler>();
      break;
    case IntegrationMethod::ImplicitEuler:
      m_timeIntegration = &LPTSolver<DEBUG_LEVEL, NDIM, P>::timeIntegration<IntegrationMethod::ImplicitEuler>;
      setForceModel<IntegrationMethod::ImplicitEuler>();
      break;
    default:
      TERMM(-1, "Invalid integration method: " + integrationName + " " + to_string(static_cast<GInt>(integrationMethod(integrationName))));
  }
}

/// Set the force model
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
template <IntegrationMethod IM>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::setForceModel() {
  GString modelName = opt_config_value("model", static_cast<GString>("constDensityRatioGravBuoStokesDrag"));

  switch(force::model(modelName)) {
    case force::Model::constDensityRatioGrav:
      m_calcA = &LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA<force::Model::constDensityRatioGrav, IM>;
      break;
    case force::Model::constDensityRatioGravStokesDrag:
      m_calcA = &LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA<force::Model::constDensityRatioGravStokesDrag, IM>;
      break;
    case force::Model::constDensityRatioGravNlinDrag:
      m_calcA = &LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA<force::Model::constDensityRatioGravNlinDrag, IM>;
      break;
    case force::Model::constDensityRatioGravBuo:
      m_calcA = &LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA<force::Model::constDensityRatioGravBuo, IM>;
      break;
    case force::Model::constDensityRatioGravBuoStokesDrag:
      m_calcA = &LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA<force::Model::constDensityRatioGravBuoStokesDrag, IM>;
      break;
    case force::Model::constDensityRatioGravBuoNlinDrag:
      m_calcA = &LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA<force::Model::constDensityRatioGravBuoNlinDrag, IM>;
      break;
    case force::Model::constDensityaGravBuo:
      m_calcA = &LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA<force::Model::constDensityaGravBuo, IM>;
      break;
    case force::Model::constDensityaGravBuoStokesDrag:
      m_calcA = &LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA<force::Model::constDensityaGravBuoStokesDrag, IM>;
      break;
    case force::Model::constDensityaGravBuoNlinDrag:
      m_calcA = &LPTSolver<DEBUG_LEVEL, NDIM, P>::calcA<force::Model::constDensityaGravBuoNlinDrag, IM>;
      break;
    default:
      TERMM(-1, "Invalid integration method: " + modelName + " " + to_string(static_cast<GInt>(integrationMethod(modelName))));
  }
}

/// Allocate the memory for the solver
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::allocateMemory() {
  // allocate memory
  m_vars.resize(m_capacity * NVARS);
  // fill variables with NAN to crash when using unset variables
  fill(m_vars.begin(), m_vars.end(), NAN);
  cerr0 << "allocated memory: " << mem() << "KB" << std::endl;
}

/// Run the main loop of the solver
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
/// \return Error code
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
auto LPTSolver<DEBUG_LEVEL, NDIM, P>::run() -> GInt {
  RECORD_TIMER_START(TimeKeeper[Timers::LPTInit]);
  loadConfiguration();
  allocateMemory();
  //  POST::init();
  initialCondition();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LPTInit]);

  RECORD_TIMER_START(TimeKeeper[Timers::LPTMainLoop]);
  // todo:implement
  //   executePostprocess(pp::HOOK::ATSTART);
  for(m_timeStep = 0; m_timeStep < m_maxNoSteps; ++m_timeStep) {
    RECORD_TIMER_START(TimeKeeper[Timers::LPTCalc]);
    timeStep();
    RECORD_TIMER_STOP(TimeKeeper[Timers::LPTCalc]);


    // also write an output if it's the last time step or if the solution has converged
    RECORD_TIMER_START(TimeKeeper[Timers::LPTIo]);
    output(m_timeStep == m_maxNoSteps - 1);
    RECORD_TIMER_STOP(TimeKeeper[Timers::LPTIo]);
  }

  if(has_config_value("analyticalSolution")) {
    compareToAnalyticalResult();
  }

  //  executePostprocess(pp::HOOK::ATEND);

  logger << "LPT Solver finished <||" << endl;
  cout << "LPT Solver finished <||" << endl;
  unusedConfigValues();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LPTMainLoop]);
  return 0;
}

/// Inital condition of the solver
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
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

/// Init particles at random position with a given volume
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::init_randomVolPos() {
  auto config = getAccessor("randomvol_pos");
  // todo: make settable
  GString       m_volType   = "box";
  GInt          noParticles = config->template required_config_value<GInt>("noparticles");
  VectorD<NDIM> cornerA     = eigenutil::unpack<NDIM>(config->template required_config_value<std::vector<GDouble>>("volume", "A"));
  VectorD<NDIM> cornerB     = eigenutil::unpack<NDIM>(config->template required_config_value<std::vector<GDouble>>("volume", "B"));

  m_init_v                  = eigenutil::unpack<NDIM>(config->template required_config_value<std::vector<GDouble>>("init_velo"));
  const GDouble initDensity = config->template required_config_value<GDouble>("init_rho_p");
  const GDouble initRadius  = config->template required_config_value<GDouble>("init_r_p");

  cerr0 << "Placing " << noParticles << " particles inside " << m_volType << " at random position" << std::endl;

  for(GInt partId = 0; partId < noParticles; ++partId) {
    center(partId)   = randomPos<NDIM>(cornerA, cornerB);
    velocity(partId) = m_init_v;
    density(partId)  = initDensity;
    radius(partId)   = initRadius;
    volume(partId)   = sphere::volumeR(initRadius);
    m_noParticles++;
  }
}

/// Perform one timestep
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::timeStep() {
  // 1. particle generation step
  generateNewParticles();
  // 2. Step update acceleration
  m_calcA(this);
  // 3. Integrate to update velocity, position
  m_timeIntegration(this);
  // 4. collision between geometry and particles
  collision();
  // 5. Delete invalid particles
  deleteInvalidParticles();
  m_currentTime += m_dt;
}

/// Calculate the acceleration of the particles
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
/// \tparam FM Force model to be used
/// \tparam IM Integration method to be used
// todo: roll all the ifs into constexpr functions
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
    switch(FM) {
      case Model::constDensityRatioGravStokesDrag:
      case Model::constDensityRatioGravBuoStokesDrag:
      case Model::constDensityaGravBuoStokesDrag:
      case Model::constDensityRatioGravNlinDrag:
      case Model::constDensityRatioGravBuoNlinDrag:
      case Model::constDensityaGravBuoNlinDrag:
        DC(partId) = 18.0 * nu_a / (4.0 * r_p * r_p * rho_p);

        // for nonlinear drag models calculate correction factor
        if(FM == Model::constDensityRatioGravNlinDrag || FM == Model::constDensityRatioGravBuoNlinDrag
           || FM == Model::constDensityaGravBuoNlinDrag) {
          re_p = rho_a * velo_p.norm() * 2.0 * r_p / nu_a;
          // todo: allow selecting drag model
          DC(partId) *= dragCoefficient<DragModel::Putnam61>(re_p);
        }

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

/// Perform the time integration
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
/// \tparam IM Integration method
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

/// Function which switches between the different modes which can generate new particles i.e. injection, breakup...
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::generateNewParticles() {
  if(m_generationMethod != GenerationMethod::None) {
    switch(m_generationMethod) {
      case GenerationMethod::ConstantRate:
        generateConst();
        break;
      case GenerationMethod::InjectionModel:
        injection();
        break;
      default:
        TERMM(-1, "Invalid generation method selected.");
    }
  }
}

/// Generate particles with a constant rate (mostly for testcases)
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
// todo: implement
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::generateConst() {}

/// Use an injection model to introduce new particles
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::injection() {
  if constexpr(NDIM == 3) {
    GInt totalNoOfParticlesInjected = 0;
    // todo: move
    // todo: make seed settable
    static randxor prng(1223);

    VectorD<NDIM> initialLocation;
    VectorD<NDIM> initialV;
    for(auto& inj : m_injectors) {
      const GInt noOfParticlesToInject = inj.timeStep(m_dt);
      for(GInt newP = 0; newP < noOfParticlesToInject; ++newP) {
        initialLocation  = inj.location() + randomPoint_inCircle(inj.orientation(), inj.holeDiameter(), prng);
        initialV         = inj.injectionVelocity() * randomNormalizedDirection_inCone(inj.orientation(), 0.5 * inj.openingAngle(), prng);
        const GInt index = addParticle(initialLocation, initialV, inj.injectedDensity(), 0.5 * inj.dropletDiameter());
      }

      totalNoOfParticlesInjected += noOfParticlesToInject;
    }
    cerr0 << "Injected " << totalNoOfParticlesInjected << std::endl;
  } else {
    TERMM(-1, "Not implemented");
  }
}

/// Add a new particle
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
/// \param pos The position of the new particle
/// \param velo The velocity of the new particle
/// \param rho The density of the new particle
/// \param r The radius of the new particle
/// \return Id of the added particle
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
auto LPTSolver<DEBUG_LEVEL, NDIM, P>::addParticle(const VectorD<NDIM>& pos, const VectorD<NDIM>& velo, const GDouble rho, const GDouble r)
    -> GInt {
  const GInt partId = m_noParticles;
  ++m_noParticles;

  if(DEBUG_LEVEL == Debug_Level::max_debug && !isnan(center(partId, 0))) {
    TERMM(-1, "overwriting");
  }

  center(partId)   = pos;
  velocity(partId) = velo;
  density(partId)  = rho;
  radius(partId)   = r;
  volume(partId)   = sphere::volumeR(r);
  return partId;
}

/// Perform the collision of particles with geometry/particles etc.
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
// todo: implement
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::collision() {}

/// Delete invalid particles e.g. particles of insufficient size, particles which have left the computational domain
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
// todo: implement
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::deleteInvalidParticles() {
  // 1. Mark particles for deletion which have left the computational domain

  // 2. Delete particles which are marked for deletion
}

/// Save the output of the particles
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
/// \param forced Perform an output *now*
/// \param postfix Add a postfix to the written files
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

/// Compare the result to an analytical solution
/// \tparam DEBUG_LEVEL Debug mode
/// \tparam NDIM Dimensionality
/// \tparam P Particle type
template <Debug_Level DEBUG_LEVEL, GInt NDIM, LPTType P>
void LPTSolver<DEBUG_LEVEL, NDIM, P>::compareToAnalyticalResult() {
  const auto analyticalSolutionName = required_config_value<GString>("analyticalSolution");

  GDouble maxError      = 0;
  GDouble sumError      = 0;
  GDouble sumErrorSq    = 0;
  GDouble sumSolution   = 0;
  GDouble sumSolutionSq = 0;

  for(GInt partId = 0; partId < m_noParticles; ++partId) {
    VectorD<NDIM> ana_sol = analytical::lpt::freefall_stokes_vel<NDIM, P>(part(partId), m_init_v, m_nu_a_infty, m_rho_a_infty,
                                                                          m_velo_a_infty, m_currentTime, m_gravity);
    const GDouble error   = (ana_sol - velocity(partId)).norm();
    maxError              = std::max(error, maxError);
    sumErrorSq += gcem::pow(error, 2);
    sumSolution += ana_sol.norm();
    sumSolutionSq += ana_sol.squaredNorm();
  }

  const GDouble L2error = (gcem::sqrt(sumErrorSq) / gcem::sqrt(sumSolutionSq)) / gcem::pow(m_noParticles, 1.0 / NDIM);
  const GDouble gre     = sumError / sumSolution;

  //  cerr0 << "Comparing to analytical result " << analyticalSolutionName << std::endl;
  //  logger << "Comparing to analytical result " << analyticalSolutionName << std::endl;
  cerr0 << "max. Error: " << maxError << std::endl;
  logger << "max. Error: " << maxError << std::endl;
  cerr0 << "avg. L2: " << L2error << std::endl;
  logger << "avg. L2: " << L2error << std::endl;
  cerr0 << "global relative error: " << gre << std::endl;
  logger << "global relative error: " << L2error << std::endl;

  // compare with the maximum expected error to give a pass/no-pass result
  const GDouble maximumExpectedError    = opt_config_value("errorMax", 1.0);
  const GDouble maximumExpectedErrorL2  = opt_config_value("errorL2", 1.0);
  const GDouble maximumExpectedErrorGRE = opt_config_value("errorGRE", 1.0);

  GBool failedErrorCheck = false;
  if(maximumExpectedError < maxError) {
    failedErrorCheck = true;
    cerr0 << "Error bounds for the maximum error have failed!!! ( < " << maximumExpectedError << ")" << std::endl;
    logger << "Error bounds for the maximum error have failed!!! ( < " << maximumExpectedError << ")" << std::endl;
  }
  if(maximumExpectedErrorL2 < L2error || isnan(L2error)) {
    failedErrorCheck = true;
    cerr0 << "Error bounds for the L2 error have failed!!! ( < " << maximumExpectedErrorL2 << ")" << std::endl;
    logger << "Error bounds for the L2 error have failed!!! ( < " << maximumExpectedErrorL2 << ")" << std::endl;
  }
  if(maximumExpectedErrorGRE < gre || isnan(gre)) {
    failedErrorCheck = true;
    cerr0 << "Error bounds for the global relative error have failed!!! ( < " << maximumExpectedErrorGRE << ")" << std::endl;
    logger << "Error bounds for the global relative error have failed!!! ( < " << maximumExpectedErrorGRE << ")" << std::endl;
  }

  if(failedErrorCheck) {
    TERMM(-1, "Analytical testcase failed");
  }
}