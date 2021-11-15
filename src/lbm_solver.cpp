#include "lbm_solver.h"
#include "pv.h"
#include "analytical_solutions.h"

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
  configuration::loadConfiguration("solver");
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMInit]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::initBenchmark(int argc, GChar** argv) {
  m_benchmark             = true;
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


  // Re = rho * u * L/nu
  // u = ma * sqrt(gamma * R * T)
  //  m_nu = m_ma / gcem::sqrt(3) / m_re * m_refLength;
  //todo: add option to define nu
  m_relaxTime = opt_config_value<GDouble>("relaxation", m_relaxTime);
  m_re        = required_config_value<GDouble>("reynoldsnumber");

  m_omega     = 1.0 / m_relaxTime;
  m_nu        = (2.0 * m_relaxTime - 1) / 6.0;
  m_refU      = m_re * m_nu / m_refLength;

  m_bndManager = std::make_unique<LBMBndManager<DEBUG_LEVEL, LBTYPE>>();

  //todo: doesn't work??
  //  std::function<const Surface<NDIM>&(GInt)> _bndrySurface = [this](const GInt id) { return bndrySurface(id); };
  using namespace placeholders;
  std::function<const Surface<NDIM>&(GString)> _bndrySurface = std::bind(&LBMSolver::bndrySurface, this, _1);
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


//template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
//void LBMSolver<DEBUG_LEVEL, LBTYPE>::setupMethod() {
//  //  switch(m_method) {
//  //    case LBMethodType::D1Q3:
//  //      m_tangentialVelo = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};
//  //      break;
//  //    case LBMethodType::D2Q5:
//  //      m_tangentialVelo = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0};
//  //      break;
//  //    case LBMethodType::D2Q9:
//  //      m_tangentialVelo = {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 4.0 / 9.0};
//  //      break;
//  //    default:
//  //      TERMM(-1, "Invalid LBM method type");
//  //  }
//}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
auto LBMSolver<DEBUG_LEVEL, LBTYPE>::run() -> GInt {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMInit]);
  loadConfiguration();
  allocateMemory();
  initialCondition();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMInit]);

  RECORD_TIMER_START(TimeKeeper[Timers::LBMMainLoop]);
  // todo: make settable
  const GInt noTimesteps = required_config_value<GInt>("maxSteps");
  GBool      converged   = false;

  for(m_timeStep = 0; m_timeStep < noTimesteps && !converged; ++m_timeStep) {

    converged = convergenceCondition();
    timeStep();

    // also write an output if its the last time step or if the solution has converged
    output( m_timeStep == noTimesteps || converged);
  }

  //todo: make settable
  static constexpr GBool analyticalTest = true;
  if(analyticalTest) {
    compareToAnalyticalResult();
  }

  logger << "LBM Solver finished <||" << endl;
  cout << "LBM Solver finished <||" << endl;
  unusedConfigValues();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMMainLoop]);
  return 0;
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
auto LBMSolver<DEBUG_LEVEL, LBTYPE>::convergenceCondition() -> GBool{
  if(m_timeStep > 0 && m_timeStep % m_outputInfoInterval == 0) {
    cerr0 <<  m_timeStep << ": ";
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

    // todo: make settable
    // todo: make independent of output
    const GDouble maxConv = std::max(std::max(std::max(convUx, convUy), convUz), convRho);
    if(m_timeStep > 1 && maxConv < 1E-12) {
      logger << "Reached convergence to: " << maxConv << std::endl;
      cerr0 << "Reached convergence to: " << maxConv << std::endl;
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
    rho<NDIM>(cellId) = 1.0;

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
  currToOldVars();
  updateMacroscopicValues();
  calcEquilibriumMoments();
  collisionStep();
  prePropBoundaryCnd();
  propagationStep();
  boundaryCnd();
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::output(const GBool forced) {

  if((m_timeStep > 0 && m_timeStep % m_outputSolutionInterval == 0) || forced) {
    std::vector<GDouble>              tmpU;
    std::vector<GDouble>              tmpV;
    std::vector<GDouble>              tmpRho;
    std::vector<GString>              index;
    std::vector<std::vector<GString>> values;

    // only output leaf cells (i.e. cells without children)
    std::function<GBool(GInt)> isLeaf = [&](GInt cellId) { return noChildren(cellId) == 0; };

    tmpU.resize(size());
    tmpV.resize(size());
    tmpRho.resize(size());
    for(GInt cellId = 0; cellId < size(); ++cellId) {
      tmpU[cellId] = velocity<NDIM>(cellId, 0);
      tmpV[cellId] = velocity<NDIM>(cellId,1);
      tmpRho[cellId] = rho<NDIM>(cellId);
    }

    index.emplace_back("u");
    index.emplace_back("v");
    index.emplace_back("rho");
    values.emplace_back(toStringVector(tmpU, size()));
    values.emplace_back(toStringVector(tmpV, size()));
    values.emplace_back(toStringVector(tmpRho, size()));

    VTK::ASCII::writePoints<NDIM>("test_" + std::to_string(m_timeStep), size(), center(), index, values, isLeaf);
  }

}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::compareToAnalyticalResult() {
  const auto analyticalSolutionName = required_config_value<GString>("analyticalSolution");
  const auto anaSolution = analytical::ns::getAnalyticalSolution(analyticalSolutionName);
  GDouble              maxError = 0;
  std::vector<GDouble> error;
  for(GInt cellId = 0; cellId < grid().size(); ++cellId) {
    const GDouble solution = anaSolution(center(cellId, 1));
    const GDouble delta    = velocity<NDIM>(cellId, 0) - solution;
    maxError               = std::max(std::abs(delta), maxError);
    error.emplace_back(std::sqrt(delta * delta / (solution * solution)));
  }
  const GDouble L2error = 1.0 / grid().size() * std::accumulate(error.begin(), error.end(), 0.0);

  const GDouble maximumExpectedError = opt_config_value("errorMax", 1.0);
  const GDouble maximumExpectedErrorL2 = opt_config_value("errorL2", 1.0);

  cerr0 << "Comparing to analytical result:" << std::endl;
  logger << "Comparing to analytical result:" << std::endl;
  cerr0 << "max. Error: " << maxError << std::endl;
  logger << "max. Error: " << maxError << std::endl;
  cerr0 << "avg. L2: " << L2error << std::endl;
  logger << "avg. L2: " << L2error << std::endl;

  GBool failedErrorCheck = false;
  if(maximumExpectedError < maxError){
    failedErrorCheck = true;
    cerr0<< "Error bounds for the maximum error have failed!!! ( < " << maximumExpectedError << ")"<< std::endl;
    logger<< "Error bounds for the maximum error have failed!!! ( < " << maximumExpectedError << ")"<< std::endl;
  }
  if(maximumExpectedErrorL2 < L2error){
    failedErrorCheck = true;
    cerr0<< "Error bounds for the L2 error have failed!!! ( < " << maximumExpectedErrorL2 << ")"<< std::endl;
    logger<< "Error bounds for the L2 error have failed!!! ( < " << maximumExpectedErrorL2 << ")"<< std::endl;
  }

  if(failedErrorCheck){
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
    for(GInt dir = 0; dir < NDIM; ++dir) {
      velocity<NDIM>(cellId, dir) = m_fold[cellId * NDIST + moment[dir][0]] + m_fold[cellId * NDIST + moment[dir][1]]
                                    + m_fold[cellId * NDIST + moment[dir][2]] - m_fold[cellId * NDIST + moment[dir][3]]
                                    - m_fold[cellId * NDIST + moment[dir][4]] - m_fold[cellId * NDIST + moment[dir][5]];
    }
    //    cerr0 << "cellId: " << cellId << " (" << velocity<NDIM>(cellId, 0) << ", " << velocity<NDIM>(cellId, 1) << ")" << std::endl;
    rho<NDIM>(cellId) = std::accumulate((m_fold.begin() + cellId * NDIST), (m_fold.begin() + (cellId + 1) * NDIST), 0.0);
  }
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::calcEquilibriumMoments() {
  static constexpr std::array<GDouble, 9> cx = {-1, 1, 0, 0, 1, 1, -1, -1, 0};
  static constexpr std::array<GDouble, 9> cy = {0, 0, -1, 1, 1, -1, -1, 1, 0};
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    for(GInt dist = 0; dist < NDIST; ++dist) {
      m_feq[cellId * NDIST + dist] =
          LBMethod<LBTYPE>::m_weights[dist]
          * (rho<NDIM>(cellId) + 3 * (velocity<NDIM>(cellId, 0) * cx[dist] + velocity<NDIM>(cellId, 1) * cy[dist]));
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

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
void LBMSolver<DEBUG_LEVEL, LBTYPE>::prePropBoundaryCnd() {
  // todo: replace by lambda function
  using namespace std::placeholders;
  std::function<GDouble&(GInt, GInt)> _fold = std::bind(&LBMSolver::fold, this, _1, _2);
  std::function<GDouble&(GInt, GInt)> _f    = std::bind(&LBMSolver::f, this, _1, _2);

  m_bndManager->preApply(_f, _fold);
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

//todo: simplify (move somewhere else?!)
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

//template class LBMSolver<Debug_Level::no_debug, LBMethodType::D1Q3>;
//template class LBMSolver<Debug_Level::min_debug, LBMethodType::D1Q3>;
//template class LBMSolver<Debug_Level::debug, LBMethodType::D1Q3>;
//template class LBMSolver<Debug_Level::more_debug, LBMethodType::D1Q3>;
//template class LBMSolver<Debug_Level::max_debug, LBMethodType::D1Q3>;
//
//template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q5>;
//template class LBMSolver<Debug_Level::min_debug, LBMethodType::D2Q5>;
//template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q5>;
//template class LBMSolver<Debug_Level::more_debug, LBMethodType::D2Q5>;
//template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q5>;

template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q9>;
template class LBMSolver<Debug_Level::min_debug, LBMethodType::D2Q9>;
template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q9>;
template class LBMSolver<Debug_Level::more_debug, LBMethodType::D2Q9>;
template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q9>;

//template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q15>;
//template class LBMSolver<Debug_Level::min_debug, LBMethodType::D3Q15>;
//template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q15>;
//template class LBMSolver<Debug_Level::more_debug, LBMethodType::D3Q15>;
//template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q15>;
