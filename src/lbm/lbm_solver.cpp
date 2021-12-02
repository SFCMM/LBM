#include "lbm_solver.h"
#include "analytical_solutions.h"
#include "lbm_pv.h"

#include <set>

using namespace std;

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::init(int argc, GChar** argv) {
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

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::init(int argc, GChar** argv, GString config_file) {
  setConfiguration(config_file);
  init(argc, argv);
  logger << "LBM Solver started ||>" << endl;
  cout << "LBM Solver started ||>" << endl;
  Configuration::load("solver");
  POST::setConfAccessor(Configuration::getAccessor("postprocessing"));
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMInit]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::initBenchmark(int argc, GChar** argv) {
  m_benchmark = true;
  init(argc, argv);

  logger << "Setting up benchmarking!" << endl;

  TERMM(-1, "Not implemented!");
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::initTimers() {
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMSolverTotal], "Total run time of the LBM Solver.", TimeKeeper[Timers::timertotal]);
  RECORD_TIMER_START(TimeKeeper[Timers::LBMSolverTotal]);

  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMInit], "Initialization of the LBM solver!", TimeKeeper[Timers::LBMSolverTotal]);
  RECORD_TIMER_START(TimeKeeper[Timers::LBMInit]);

  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMMainLoop], "Main Loop of the LBM solver!", TimeKeeper[Timers::LBMSolverTotal]);
}


template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::loadConfiguration() {
  m_solverType = LBSolverType::BGK; // todo: load from file

  m_outputInfoInterval     = opt_config_value<GInt>("info_interval", m_outputInfoInterval);
  m_outputSolutionInterval = opt_config_value<GInt>("solution_interval", m_outputSolutionInterval);


  //  const GDouble m_referenceLength = 1.0;
  //  const GDouble m_ma = 0.001;

  // Re = rho * u * L/nu
  // u = ma * sqrt(gamma * R * T)
  //  m_nu = m_ma / gcem::sqrt(3) / m_re * m_refLength;
  // todo: add option to define nu
  m_relaxTime = opt_config_value<GDouble>("relaxation", m_relaxTime);
  m_re        = required_config_value<GDouble>("reynoldsnumber");


  /// todo: 1.73205080756887729352??? (F1BCS)
  //  m_nu    = m_ma / 1.73205080756887729352 / m_re * m_referenceLength; //* (FFPOW2[maxLevel() - a_level(pCellId)]);
  //  m_omega = 2.0 / (1.0 + 6.0 * m_nu);

  m_omega = 1.0 / m_relaxTime;
  m_nu    = (2.0 * m_relaxTime - 1) / 6.0;
  m_refU  = m_re * m_nu / m_refLength;

  m_bndManager = std::make_unique<LBMBndManager<DEBUG_LEVEL, LBTYPE>>();

  // todo: doesn't work??
  //   std::function<const Surface<NDIM>&(GInt)> _bndrySurface = [this](const GInt id) { return bndrySurface(id); };
  using namespace placeholders;
  std::function<const Surface<DEBUG_LEVEL, NDIM>&(GString)> _bndrySurface = std::bind(&LBMSolver::bndrySurface, this, _1);
  m_bndManager->setupBndryCnds(config()["boundary"], _bndrySurface);


  cerr0 << "<<<<<<<<<<<<>>>>>>>>>>>>>" << std::endl;
  cerr0 << "LBM Type "
        << "BGK" << std::endl; // todo: fix me
  cerr0 << "LBM Method "
        << "D2Q9" << std::endl; // todo: fix me
  cerr0 << "No. of variables " << std::to_string(NVARS) << std::endl;
  cerr0 << "No. Leaf cells: " << noLeafCells() << std::endl;
  cerr0 << "No. Bnd cells: " << noBndCells() << std::endl;
  cerr0 << "Relaxation Time: " << m_relaxTime << std::endl;
  cerr0 << "Omega: " << m_omega << std::endl;
  cerr0 << "Reynolds Number: " << m_re << std::endl;
  cerr0 << "Viscosity: " << m_nu << std::endl;
  cerr0 << "Reference V: " << m_refU << std::endl;
  cerr0 << "+++++++++++++++++++++++++" << std::endl;
}


template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::allocateMemory() {
  // allocate memory
  m_f.resize(grid().size() * NDIST);
  m_feq.resize(grid().size() * NDIST);
  m_fold.resize(grid().size() * NDIST);
  m_vars.resize(grid().size() * NVARS);
  m_varsold.resize(grid().size() * NVARS);
}


// template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
// void LBMSolver<DEBUG_LEVEL, LBTYPE>::setupMethod() {
//   //  switch(m_method) {
//   //    case LBMethodType::D1Q3:
//   //      m_tangentialVelo = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};
//   //      break;
//   //    case LBMethodType::D2Q5:
//   //      m_tangentialVelo = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0};
//   //      break;
//   //    case LBMethodType::D2Q9:
//   //      m_tangentialVelo = {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 4.0 / 9.0};
//   //      break;
//   //    default:
//   //      TERMM(-1, "Invalid LBM method type");
//   //  }
// }

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
auto LBMSolver<DEBUG_LEVEL, LBTYPE>::run() -> GInt {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMInit]);
  loadConfiguration();
  allocateMemory();
  POST::init();
  initialCondition();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMInit]);

  RECORD_TIMER_START(TimeKeeper[Timers::LBMMainLoop]);
  // todo: make settable
  const GInt noTimesteps = required_config_value<GInt>("maxSteps");
  GBool      converged   = false;

  executePostprocess(pp::HOOK::ATSTART);

  for(m_timeStep = 0; m_timeStep < noTimesteps && !converged; ++m_timeStep) {
    converged = convergenceCondition();
    timeStep();

    // also write an output if its the last time step or if the solution has converged
    output(m_timeStep == noTimesteps || converged);
  }

  // todo: make settable
  static constexpr GBool analyticalTest = true;
  if(analyticalTest) {
    compareToAnalyticalResult();
  }

  executePostprocess(pp::HOOK::ATEND);

  logger << "LBM Solver finished <||" << endl;
  cout << "LBM Solver finished <||" << endl;
  unusedConfigValues();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMMainLoop]);
  return 0;
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
auto LBMSolver<DEBUG_LEVEL, LBTYPE>::convergenceCondition() -> GBool {
  if(m_timeStep > 0 && m_timeStep % m_outputInfoInterval == 0) {
    cerr0 << m_timeStep << ": ";
    logger << m_timeStep << ": ";

    const GDouble convUx  = sumAbsDiff(PV::velocitiy(0));
    const GDouble convUy  = sumAbsDiff(PV::velocitiy(1));
    const GDouble convUz  = (NDIM == 3) ? sumAbsDiff(PV::velocitiy(2)) : NAN;
    const GDouble convRho = sumAbsDiff(PV::rho<NDIM>());

    cerr0 << "u_x: " << convUx << " ";
    cerr0 << "u_y: " << convUy << " ";
    if(NDIM == 3) {
      cerr0 << "u_z: " << convUz << " ";
    }
    cerr0 << "rho: " << convRho << " ";
    cerr0 << std::endl;

    logger << "u_x: " << convUx << " ";
    logger << "u_y: " << convUy << " ";
    if(NDIM == 3) {
      logger << "u_z: " << convUz << " ";
    }
    logger << "rho: " << convRho << " ";
    logger << std::endl;

    // todo: make independent of output
    const GDouble maxConv         = std::max(std::max(std::max(convUx, convUy), convUz), convRho);
    const GDouble convergenceCrit = opt_config_value("convergence", 1E-12);
    if(m_timeStep > 1 && maxConv < convergenceCrit) {
      logger << "Reached convergence to: " << maxConv << std::endl;
      cerr0 << "Reached convergence to: " << maxConv << std::endl;
      return true;
    }
    if(m_timeStep > 1 && std::isnan(maxConv)) {
      logger << "Solution diverged!" << std::endl;
      cerr0 << "Solution diverged!" << std::endl;
      m_diverged = true;
      // return true to cancel further calculation
      return true;
    }
  }
  return false;
}


template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::initialCondition() {
  // init to zero:
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    for(const auto dir : PV::velocities<NDIM>()) {
      m_vars[cellId * NVARS + dir]    = 0;
      m_varsold[cellId * NVARS + dir] = 0;
    }
    rho(cellId) = 1.0;

    // assuming initial zero velocity and density 1
    for(GInt dist = 0; dist < NDIST; ++dist) {
      m_feq[cellId * NDIST + dist]  = LBMethod<LBTYPE>::m_weights[dist];
      m_f[cellId * NDIST + dist]    = m_feq[cellId * NDIST + dist];
      m_fold[cellId * NDIST + dist] = m_feq[cellId * NDIST + dist];
    }
  }
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::timeStep() {
  executePostprocess(pp::HOOK::BEFORETIMESTEP);
  currToOldVars();
  updateMacroscopicValues();
  calcEquilibriumMoments();
  collisionStep();
  forcing();
  prePropBoundaryCnd();
  propagationStep();
  boundaryCnd();
  executePostprocess(pp::HOOK::AFTERTIMESTEP);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::output(const GBool forced, const GString& postfix) {
  if(m_diverged) {
    // output the values before updating the macroscopic values in case the solution diverged
    output(true, "bdiv");
  }
  //  updateMacroscopicValues();
  if((m_timeStep > 0 && m_timeStep % m_outputSolutionInterval == 0) || forced) {
    std::array<std::vector<GDouble>, NDIM> tmpVel;
    std::vector<GDouble>                   tmpRho;
    std::vector<IOIndex>                   index;
    std::vector<std::vector<GString>>      values;

    // only output leaf cells (i.e. cells without children)
    std::function<GBool(GInt)> isLeaf = [&](GInt cellId) { return noChildren(cellId) == 0; };

    for(GInt dir = 0; dir < NDIM; ++dir) {
      tmpVel[dir].resize(size());
    }
    tmpRho.resize(size());

    for(GInt cellId = 0; cellId < size(); ++cellId) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        tmpVel[dir][cellId] = velocity(cellId, dir);
      }
      tmpRho[cellId] = rho(cellId);
    }

    for(GInt dir = 0; dir < NDIM; ++dir) {
      // todo: fix type
      index.emplace_back(IOIndex{static_cast<GString>(PV::VELSTR[dir]), "float32"});
      values.emplace_back(toStringVector(tmpVel[dir], size()));
    }
    // todo: fix type
    index.emplace_back(IOIndex{"rho", "float32"});
    values.emplace_back(toStringVector(tmpRho, size()));

    VTK::ASCII::writePoints<NDIM>("test_" + std::to_string(m_timeStep) + postfix, size(), center(), index, values, isLeaf);
  }
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::compareToAnalyticalResult() {
  // for an analytical testcase the value "analyticalSolution" needs to be defined
  const auto analyticalSolutionName = required_config_value<GString>("analyticalSolution");
  // from the analyticalSolutionName we can get a functional representation of the solution
  const auto anaSolution = analytical::ns::getAnalyticalSolution(analyticalSolutionName);

  // determine the error
  GDouble              maxError = 0;
  std::vector<GDouble> error;

  // exclude cells from error calculation, for example, the boundary surfaces
  std::vector<std::vector<GInt>> m_excludedCells;
  if(has_config_value("analyticalSolutionExcludeSurface")) {
    const auto excludedSurfaces = required_config_value<std::vector<GString>>("analyticalSolutionExcludeSurface");
    for(const auto& exclSurfName : excludedSurfaces) {
      m_excludedCells.emplace_back(bndrySurface(exclSurfName).getCellList());
    }
  }
  for(GInt cellId = 0; cellId < grid().size(); ++cellId) {
    GBool excluded = false;
    if(!m_excludedCells.empty()) {
      for(const auto& cellList : m_excludedCells) {
        if(std::find(cellList.begin(), cellList.end(), cellId) != cellList.end()) {
          excluded = true;
          break;
        }
      }
    }
    if(excluded) {
      continue;
    }
    const GDouble solution = anaSolution(center(cellId, 1));
    const GDouble delta    = velocity(cellId, 0) - solution;
    maxError               = std::max(std::abs(delta), maxError);
    error.emplace_back(std::sqrt(delta * delta / (solution * solution)));
  }
  const GDouble L2error = 1.0 / grid().size() * std::accumulate(error.begin(), error.end(), 0.0);

  // if that compare with the maximum expected error to give a pass/no-pass result
  const GDouble maximumExpectedError   = opt_config_value("errorMax", 1.0);
  const GDouble maximumExpectedErrorL2 = opt_config_value("errorL2", 1.0);

  cerr0 << "Comparing to analytical result:" << std::endl;
  logger << "Comparing to analytical result:" << std::endl;
  cerr0 << "max. Error: " << maxError << std::endl;
  logger << "max. Error: " << maxError << std::endl;
  cerr0 << "avg. L2: " << L2error << std::endl;
  logger << "avg. L2: " << L2error << std::endl;

  GBool failedErrorCheck = false;
  if(maximumExpectedError < maxError) {
    failedErrorCheck = true;
    cerr0 << "Error bounds for the maximum error have failed!!! ( < " << maximumExpectedError << ")" << std::endl;
    logger << "Error bounds for the maximum error have failed!!! ( < " << maximumExpectedError << ")" << std::endl;
  }
  if(maximumExpectedErrorL2 < L2error) {
    failedErrorCheck = true;
    cerr0 << "Error bounds for the L2 error have failed!!! ( < " << maximumExpectedErrorL2 << ")" << std::endl;
    logger << "Error bounds for the L2 error have failed!!! ( < " << maximumExpectedErrorL2 << ")" << std::endl;
  }

  if(failedErrorCheck) {
    TERMM(-1, "Analytical testcase failed");
  }
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::currToOldVars() {
  std::copy(m_vars.begin(), m_vars.end(), m_varsold.begin());
  //  std::copy(m_f.begin(), m_f.end(), m_fold.begin());
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::updateMacroscopicValues() {
  static constexpr std::array<std::array<GInt, 6>, 2> moment = {{{{1, 4, 5, 0, 6, 7}}, {{3, 4, 7, 2, 5, 6}}}};

  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    // todo:skip non-leaf cells!
    rho(cellId) = std::accumulate((m_fold.begin() + cellId * NDIST), (m_fold.begin() + (cellId + 1) * NDIST), 0.0);

    for(GInt dir = 0; dir < NDIM; ++dir) {
      velocity(cellId, dir) =
          (m_fold[cellId * NDIST + moment[dir][0]] + m_fold[cellId * NDIST + moment[dir][1]] + m_fold[cellId * NDIST + moment[dir][2]]
           - m_fold[cellId * NDIST + moment[dir][3]] - m_fold[cellId * NDIST + moment[dir][4]] - m_fold[cellId * NDIST + moment[dir][5]])
          / rho(cellId);
    }
  }
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::calcEquilibriumMoments() {
  static constexpr std::array<GDouble, 9> cx = {-1, 1, 0, 0, 1, 1, -1, -1, 0};
  static constexpr std::array<GDouble, 9> cy = {0, 0, -1, 1, 1, -1, -1, 1, 0};
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    //    const GDouble squaredV = gcem::pow(velocity<NDIM>(cellId, 0), 2) + gcem::pow(velocity<NDIM>(cellId, 1), 2);
    for(GInt dist = 0; dist < NDIST; ++dist) {
      m_feq[cellId * NDIST + dist] =
          LBMethod<LBTYPE>::m_weights[dist] * rho(cellId) * (1 + 3 * (velocity(cellId, 0) * cx[dist] + velocity(cellId, 1) * cy[dist]));
      //      const GDouble cu             = velocity<NDIM>(cellId, 0) * cx[dist] + velocity<NDIM>(cellId, 1) * cy[dist];
      //      m_feq[cellId * NDIST + dist] = LBMethod<LBTYPE>::m_weights[dist] * rho<NDIM>(cellId) * (1 + 3 * cu + 4.5 * cu * cu - 1.5 *
      //      squaredV);
    }
  }
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::collisionStep() {
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    for(GInt dist = 0; dist < NDIST; ++dist) {
      m_f[cellId * NDIST + dist] = (1 - m_omega) * m_fold[cellId * NDIST + dist] + m_omega * m_feq[cellId * NDIST + dist];
    }
  }
}

// todo:refactor
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::forcing() {
  static GBool  info        = true;
  const GString forcingType = opt_config_value<GString>("forcing", "");

  if(!forcingType.empty()) {
    if(info) {
      cerr0 << "Using forcing!" << std::endl;
      logger << "Using forcing!" << std::endl;
      info = false;
    }

    const auto inlet  = bndrySurface("-x");
    const auto outlet = bndrySurface("+x");

    // const GDouble umax = 0.1;
    //    const GDouble gradP = 0.00011276015; //gradP=8*nu*u_max/(NY)^2;
    const GDouble outletPressure = 1.0;
    //    const GDouble inletPressure  = 1.010487; //rho_inlet=3*(NX-1)*gradP+outletPressure; //0.092
    //    const GDouble inletPressure  = 1.011497487;//0.09927
    const GDouble inletPressure = 1.0115985357; // 0.01


    // set Outlet forcing
    for(const GInt inletCellId : inlet.getCellList()) {
      const GInt           valueCellId = grid().neighbor(inletCellId, 1);
      const VectorD<NDIM>& centerInlet = grid().center(valueCellId);
      for(const GInt outletCellId : outlet.getCellList()) {
        const VectorD<NDIM>& centerOutlet = grid().center(outletCellId);
        if(std::abs(centerInlet[1] - centerOutlet[1]) < GDoubleEps) {
          for(GInt dist = 0; dist < LBMethod<LBTYPE>::m_noDists; ++dist) {
            f(outletCellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * outletPressure
                                        * (1
                                           + 3
                                                 * (vars(valueCellId, PV::velocitiy(0)) * LBMethod<LBTYPE>::m_dirs[dist][0]
                                                    + 0 * LBMethod<LBTYPE>::m_dirs[dist][1]))
                                    + f(valueCellId, dist) - feq(valueCellId, dist);
          }
        }
      }
    }

    // set Inlet forcing
    for(const GInt outletCellId : outlet.getCellList()) {
      const GInt           valueCellId  = grid().neighbor(outletCellId, 0);
      const VectorD<NDIM>& centerOutlet = grid().center(valueCellId);
      for(const GInt inletCellId : inlet.getCellList()) {
        const VectorD<NDIM>& centerInlet = grid().center(inletCellId);
        if(std::abs(centerInlet[1] - centerOutlet[1]) < GDoubleEps) {
          for(GInt dist = 0; dist < LBMethod<LBTYPE>::m_noDists; ++dist) {
            f(inletCellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * inletPressure
                                       * (1
                                          + 3
                                                * (vars(valueCellId, PV::velocitiy(0)) * LBMethod<LBTYPE>::m_dirs[dist][0]
                                                   + 0 * LBMethod<LBTYPE>::m_dirs[dist][1]))
                                   + f(valueCellId, dist) - feq(valueCellId, dist);
          }
        }
      }
    }
  }
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::prePropBoundaryCnd() {
  // todo: replace by lambda function
  using namespace std::placeholders;
  std::function<GDouble&(GInt, GInt)> _fold = std::bind(&LBMSolver::fold, this, _1, _2);
  std::function<GDouble&(GInt, GInt)> _f    = std::bind(&LBMSolver::f, this, _1, _2);
  std::function<GDouble&(GInt, GInt)> _v    = std::bind(&LBMSolver::vars, this, _1, _2);

  m_bndManager->preApply(_f, _fold, _v);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::propagationStep() {
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    for(GInt dist = 0; dist < NDIST - 1; ++dist) {
      const GInt neighborId = m_grid->neighbor(cellId, dist);
      if(neighborId != INVALID_CELLID) {
        m_fold[neighborId * NDIST + dist] = m_f[cellId * NDIST + dist];
      }
    }
    fold(cellId, NDIST - 1) = f(cellId, NDIST - 1);
  }
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::boundaryCnd() {
  // todo: replace by lambda function
  using namespace std::placeholders;
  std::function<GDouble&(GInt, GInt)> _fold = std::bind(&LBMSolver::fold, this, _1, _2);
  std::function<GDouble&(GInt, GInt)> _f    = std::bind(&LBMSolver::f, this, _1, _2);
  std::function<GDouble&(GInt, GInt)> _v    = std::bind(&LBMSolver::vars, this, _1, _2);

  m_bndManager->apply(_f, _fold, _v);
}

// todo: simplify (move somewhere else?!)
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::transferGrid(const GridInterface& grid) {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMInit]);
  cerr0 << "Transferring " + std::to_string(NDIM) + "D Grid to LBM solver" << std::endl;
  logger << "Transferring " + std::to_string(NDIM) + "D Grid to LBM solver" << std::endl;


  switch(NDIM) {
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

  switch(NDIM) {
    case 1:
      static_cast<CartesianGrid<DEBUG_LEVEL, 1>*>(m_grid.get())
          ->loadGridInplace(*static_cast<const CartesianGridGen<DEBUG_LEVEL, 1>*>(static_cast<const void*>(&grid)), config());
      break;
    case 2:
      static_cast<CartesianGrid<DEBUG_LEVEL, 2>*>(m_grid.get())
          ->loadGridInplace(*static_cast<const CartesianGridGen<DEBUG_LEVEL, 2>*>(static_cast<const void*>(&grid)), config());
      break;
    case 3:
      static_cast<CartesianGrid<DEBUG_LEVEL, 3>*>(m_grid.get())
          ->loadGridInplace(*static_cast<const CartesianGridGen<DEBUG_LEVEL, 3>*>(static_cast<const void*>(&grid)), config());
      break;
    default:
      TERMM(-1, "Only dimensions 1,2 and 3 are supported.");
  }

  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMInit]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
constexpr auto LBMSolver<DEBUG_LEVEL, LBTYPE>::isThermal() -> GBool {
  return false;
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
auto LBMSolver<DEBUG_LEVEL, LBTYPE>::sumAbsDiff(const GInt var) const -> GDouble {
  GDouble conv = 0.0;
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    conv += std::abs(m_vars[cellId * NVARS + var] - m_varsold[cellId * NVARS + var]);
  }
  return conv;
}

// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D1Q3>;
// template class LBMSolver<Debug_Level::min_debug, LBMethodType::D1Q3>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D1Q3>;
// template class LBMSolver<Debug_Level::more_debug, LBMethodType::D1Q3>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D1Q3>;
//
// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q5>;
// template class LBMSolver<Debug_Level::min_debug, LBMethodType::D2Q5>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q5>;
// template class LBMSolver<Debug_Level::more_debug, LBMethodType::D2Q5>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q5>;

template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q9>;
template class LBMSolver<Debug_Level::min_debug, LBMethodType::D2Q9>;
template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q9>;
template class LBMSolver<Debug_Level::more_debug, LBMethodType::D2Q9>;
template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q9>;

// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q15>;
// template class LBMSolver<Debug_Level::min_debug, LBMethodType::D3Q15>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q15>;
// template class LBMSolver<Debug_Level::more_debug, LBMethodType::D3Q15>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q15>;
