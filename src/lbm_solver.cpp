#include "lbm_solver.h"
#include "pv.h"

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
  m_noDists    = LBMethod::noDists.at(static_cast<GInt>(m_method));
  m_noVars     = m_dim + 1 + static_cast<GInt>(isThermal());
}


template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::finishInit() {
  // allocate memory
  m_f.resize(grid().size() * m_noDists);
  m_feq.resize(grid().size() * m_noDists);
  m_fold.resize(grid().size() * m_noDists);
  m_vars.resize(grid().size() * m_noVars);
  m_varsold.resize(grid().size() * m_noVars);
  setupMethod();
}


template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::setupMethod() {
  switch(m_method) {
    case LBMethodType::D1Q3:
      m_weight = {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0};
      break;
    case LBMethodType::D2Q5:
      m_weight = {1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
      break;
    case LBMethodType::D2Q9:
      m_weight = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
      break;
    default:
      TERMM(-1, "Invalid LBM method type");
  }
}


template <Debug_Level DEBUG_LEVEL>
auto LBMSolver<DEBUG_LEVEL>::run() -> GInt {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMMainLoop]);
  switch(m_dim) {
    case 1:
      initialCondition<1>();
      break;
    case 2:
      initialCondition<2>();
      break;
    case 3:
      initialCondition<3>();
      break;
    case 4:
      initialCondition<4>();
      break;
    default:
      TERMM(-1, "Invalid number of dimensions 1-4.");
  }
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
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::initialCondition() {
  // init to zero:
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    for(const auto dir : PV::velocities<NDIM>()) {
      m_vars[cellId * m_noVars + dir] = 0;
    }
    rho<NDIM>(cellId) = 1.0;

    // assuming initial zero velocity and density 1
    for(GInt dist = 0; dist < m_noDists; dist++) {
      m_feq[cellId * m_noDists + dist]  = m_weight[dist];
      m_f[cellId * m_noDists + dist]    = m_feq[cellId * m_noDists + dist];
      m_fold[cellId * m_noDists + dist] = m_feq[cellId * m_noDists + dist];
    }
  }

  // Re = rho * u * L/nu
  // u = ma * sqrt(gamma * R * T)
  //  m_nu = m_ma / gcem::sqrt(3) / m_re * m_refLength;
  // todo: make settable
  m_relaxTime = 0.9;
  m_omega = 1.0/m_relaxTime;
  m_re        = 10;
  m_nu        = (2.0 * m_relaxTime - 1) / 6.0; // default=0.133
  m_refU      = m_re * m_nu / m_refLength;     // default=1.3333
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::timeStep() {
  currToOldVars<NDIM>();
  updateMacroscopicValues<NDIM>();
  calcEquilibriumMoments<NDIM>();
  collisionStep<NDIM>();
  //  m_lbm->collisionStep();
  //  boundaryCnd<NDIM>();
  propagationStep<NDIM>();
  boundaryCnd<NDIM>();
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::currToOldVars() {
  std::copy(m_vars.begin(), m_vars.end(), m_varsold.begin());
  std::copy(m_f.begin(), m_f.end(), m_fold.end());
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::updateMacroscopicValues() {
  static constexpr std::array<std::array<GInt, 6>, 2> moment = {{{{1, 4, 5, 0, 6, 7}}, {{3, 4, 7, 2, 5, 6}}}};

  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    // todo:skip non-leaf cells!
    for(GInt dir = 0; dir < NDIM; ++dir) {
      velocity<NDIM>(cellId, dir) = m_fold[cellId * m_noDists + moment[dir][0]] + m_fold[cellId * m_noDists + moment[dir][1]]
                                    + m_fold[cellId * m_noDists + moment[dir][2]] - m_fold[cellId * m_noDists + moment[dir][3]]
                                    - m_fold[cellId * m_noDists + moment[dir][4]] - m_fold[cellId * m_noDists + moment[dir][5]];
    }
    rho<NDIM>(cellId) = std::accumulate((m_fold.begin() + cellId * m_noDists), (m_fold.begin() + (cellId + 1) * m_noDists), 0.0);
  }
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::calcEquilibriumMoments() {
  static constexpr std::array<GDouble, 9> cx = {-1, 1, 0, 0, 1, 1, -1, -1, 0};
  static constexpr std::array<GDouble, 9> cy = {0, 0, -1, 1, 1, -1, -1, 1, 0};
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    for(GInt dist = 0; dist < m_noDists; ++dist) {
      m_feq[cellId * m_noDists + dist] =
          m_weight[dist] * (rho<NDIM>(cellId) + 3 * (velocity<NDIM>(cellId, 0) * cx[dist] + velocity<NDIM>(cellId, 1) * cy[dist]));
    }
  }
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::collisionStep() {
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    for(GInt dist = 0; dist < m_noDists; ++dist) {
      m_f[cellId * m_noDists + dist] = (1-m_omega) * m_fold[cellId * m_noDists + dist] + m_omega * m_feq[cellId * m_noDists + dist];
    }
  }
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::propagationStep() {
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    for(GInt dist = 0; dist < m_noDists; ++dist) {
      const GInt neighborId = m_grid->neighbor(cellId, dist);
      if(neighborId != INVALID_CELLID) {
        m_fold[neighborId] = m_f[cellId];
      }
    }
  }
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::boundaryCnd() {}

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
constexpr auto LBMSolver<DEBUG_LEVEL>::isThermal() -> GBool {
  return false;
}


template class LBMSolver<Debug_Level::no_debug>;
template class LBMSolver<Debug_Level::min_debug>;
template class LBMSolver<Debug_Level::debug>;
template class LBMSolver<Debug_Level::more_debug>;
template class LBMSolver<Debug_Level::max_debug>;