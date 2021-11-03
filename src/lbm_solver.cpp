#include "lbm_solver.h"
#include "pv.h"

#include <set>

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
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::loadConfiguration() {
  m_solverType         = LBSolverType::BGK;  // todo: load from file
  m_method             = LBMethodType::D2Q9; // todo: load from file
  m_noDists            = noDists(m_method);
  m_noVars             = m_dim + 1 + static_cast<GInt>(isThermal());
  // todo: make settable
  m_outputInfoInterval = 10;
  m_outputSolutionInterval = 100;


  // Re = rho * u * L/nu
  // u = ma * sqrt(gamma * R * T)
  //  m_nu = m_ma / gcem::sqrt(3) / m_re * m_refLength;
  // todo: make settable
  m_relaxTime = 0.9;
  m_omega     = 1.0 / m_relaxTime;
  m_re        = 0.75;
  m_nu        = (2.0 * m_relaxTime - 1) / 6.0; // default=0.133
  m_refU      = m_re * m_nu / m_refLength;     // default=1.3333

  m_bndManager = std::make_unique<LBMBndManager<DEBUG_LEVEL>>(m_dim, m_noDists);

  //todo: just for current testcase
  m_bndManager->template addBndry<NDIM>(BndryType::Wall_BounceBack, grid<NDIM>()->bndrySurface(static_cast<GInt>(LBMDir::mY)));
  m_bndManager->template addBndry<NDIM>(BndryType::Wall_BounceBack_TangentialVelocity, grid<NDIM>()->bndrySurface(static_cast<GInt>(LBMDir::pY)));



  cerr0 << "<<<<<<<<<<<<>>>>>>>>>>>>>" << std::endl;
  cerr0 << "LBM Type "
        << "BGK" << std::endl;//todo: fix me
  cerr0 << "LBM Method "
        << "D2Q9" << std::endl;//todo: fix me
  cerr0 << "No. of variables " << std::to_string(m_noVars) << std::endl;
  cerr0 << "No. Leaf cells: " << grid<NDIM>()->noLeafCells() << std::endl;
  cerr0 << "No. Bnd cells: " << grid<NDIM>()->noBndCells() << std::endl;
  cerr0 << "Relaxation Time: " << m_relaxTime << std::endl;
  cerr0 << "Reynolds Number: " << m_re << std::endl;
  cerr0 << "Viscosity: " << m_nu << std::endl;
  cerr0 << "Reference V: " << m_refU << std::endl;
  cerr0 << "+++++++++++++++++++++++++" << std::endl;
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
      m_weight = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};
      break;
    case LBMethodType::D2Q5:
      m_weight = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0};
      break;
    case LBMethodType::D2Q9:
      m_weight = {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 4.0 / 9.0};
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
      return run<1>();
      break;
    case 2:
      return run<2>();
      break;
    case 3:
      return run<3>();
      break;
    case 4:
      return run<4>();
      break;
    default:
      TERMM(-1, "Invalid number of dimensions 1-4.");
  }
  return -1;
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
auto LBMSolver<DEBUG_LEVEL>::run() -> GInt {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMMainLoop]);
  loadConfiguration<NDIM>();
  finishInit();
  initialCondition<NDIM>();
  setupPeriodicConnection<NDIM>();

  // todo: make settable
  const GInt noTimesteps = 20000;
  GBool converged = false;
  for(GInt ts = 0; ts < noTimesteps && !converged; ++ts) {
    if(m_timeStep > 0 && m_timeStep % m_outputInfoInterval == 0) {
      cerr0 << "Running LBM at " << m_timeStep << std::endl;
      logger << "Running LBM at " << m_timeStep << std::endl;

      const GDouble convUx  = convergence<NDIM>(PV::velocitiy(0));
      const GDouble convUy  = convergence<NDIM>(PV::velocitiy(1));
      const GDouble convUz  = (NDIM == 3) ? convergence<NDIM>(PV::velocitiy(2)) : NAN;
      const GDouble convRho = convergence<NDIM>(PV::rho<NDIM>());

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
      const GDouble maxConv = std::max(std::max(std::max(convUx, convUy), convUz), convRho);
      if(m_timeStep > 1 && maxConv < 1E-12){
        logger << "Reached convergence to: " << maxConv <<std::endl;
        cerr0 << "Reached convergence to: " << maxConv <<std::endl;
        converged = true;
      }
    }
    timeStep<NDIM>();

    std::vector<GDouble>              tmp;
    std::vector<GString>              index;
    std::vector<std::vector<GString>> values;

    if((m_timeStep > 0 && m_timeStep%m_outputSolutionInterval==0) || m_timeStep == noTimesteps || converged) {
      // only output leaf cells (i.e. cells without children)
      std::function<GBool(GInt)> isLeaf = [&](GInt cellId) { return grid<NDIM>()->noChildren(cellId) == 0; };

      tmp.resize(grid<NDIM>()->size());
      for(GInt cellId = 0; cellId < grid<NDIM>()->size(); ++cellId) {
        tmp[cellId] = velocity<NDIM>(cellId, 0);
      }

      index.emplace_back("u");
      values.emplace_back(toStringVector(tmp, grid<NDIM>()->size()));

      VTK::ASCII::writePoints<NDIM>("test_" + std::to_string(m_timeStep), grid<NDIM>()->size(), grid<NDIM>()->center(), index, values,
                                    isLeaf);
    }

    ++m_timeStep;
  }

  const GDouble couette_channelHeight = 5.0;
  auto analytical = [&](const GDouble y){
    return m_refU/couette_channelHeight * y;
  };

  GDouble maxError = 0;
  std::vector<GDouble> error;
  for(GInt cellId = 0; cellId < grid().size(); ++cellId){
    const GDouble solution = analytical(grid<NDIM>()->center(cellId, 1));
    const GDouble delta = velocity<NDIM>(cellId, 0) - solution;
    maxError = std::max(std::abs(delta), maxError);
    error.emplace_back(std::sqrt(delta*delta/(solution * solution)));
  }
  const GDouble L2error = 1.0/grid().size() * std::accumulate(error.begin(), error.end(), 0.0);

  cerr0 << "max. Error: " << maxError <<std::endl;
  cerr0 << "avg. L2: " << L2error <<std::endl;

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
      m_vars[cellId * m_noVars + dir]    = 0;
      m_varsold[cellId * m_noVars + dir] = 0;
    }
    rho<NDIM>(cellId) = 1.0;

    // assuming initial zero velocity and density 1
    for(GInt dist = 0; dist < m_noDists; ++dist) {
      m_feq[cellId * m_noDists + dist]  = m_weight.at(dist);
      m_f[cellId * m_noDists + dist]    = m_feq[cellId * m_noDists + dist];
      m_fold[cellId * m_noDists + dist] = m_feq[cellId * m_noDists + dist];
    }
  }
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
  //  std::copy(m_f.begin(), m_f.end(), m_fold.begin());
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
//    cerr0 << "cellId: " << cellId << " (" << velocity<NDIM>(cellId, 0) << ", " << velocity<NDIM>(cellId, 1) << ")" << std::endl;
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
      m_f[cellId * m_noDists + dist] = (1 - m_omega) * m_fold[cellId * m_noDists + dist] + m_omega * m_feq[cellId * m_noDists + dist];
      //      cerr0<< m_f[cellId * m_noDists + dist] << " m_omega " << m_omega << " m_fold " << m_fold[cellId * m_noDists + dist] << " m_feq
      //      " <<
      //          m_feq[cellId * m_noDists + dist] <<std::endl;
    }
  }
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::propagationStep() {
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    for(GInt dist = 0; dist < m_noDists - 1; ++dist) {
      const GInt neighborId = m_grid->neighbor(cellId, dist);
      if(neighborId != INVALID_CELLID) {
        m_fold[neighborId * m_noDists + dist] = m_f[cellId * m_noDists + dist];
      }
    }
  }
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::boundaryCnd() {
    using namespace  std::placeholders;
    std::function<GDouble&(GInt, GInt)> _fold = std::bind(&LBMSolver::fold, this, _1, _2);
    std::function<GDouble&(GInt, GInt)> _f = std::bind(&LBMSolver::f, this, _1, _2);

    m_bndManager->apply(_f, _fold);


//  //  cerr0 << "bnd" << std::endl;
//  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
//    //    cerr0 << "bnd " << grid().property(cellId, CellProperties::bndry) << std::endl;
//    //    cerr0 << "leaf " << grid().property(cellId, CellProperties::leaf) << std::endl;
//
//
//
//    if(grid().property(cellId, CellProperties::bndry) && grid().property(cellId, CellProperties::leaf)) {
//      //      // -X
//      //      if(grid().neighbor(cellId, 0) == INVALID_CELLID) {
//      //        m_fold[cellId * m_noDists + 1] = m_f[cellId * m_noDists + 0];
//      //        m_fold[cellId * m_noDists + 4] = m_f[cellId * m_noDists + 6];
//      //        m_fold[cellId * m_noDists + 5] = m_f[cellId * m_noDists + 7];
//      //      }
//      //      // +X
//      //      if(grid().neighbor(cellId, 1) == INVALID_CELLID) {
//      //        m_fold[cellId * m_noDists + 0] = m_f[cellId * m_noDists + 1];
//      //        m_fold[cellId * m_noDists + 6] = m_f[cellId * m_noDists + 4];
//      //        m_fold[cellId * m_noDists + 7] = m_f[cellId * m_noDists + 5];
//      //      }
//      // -y
////      if(grid().neighbor(cellId, 2) == INVALID_CELLID) {
////        m_fold[cellId * m_noDists + 3] = m_f[cellId * m_noDists + 2];
////        m_fold[cellId * m_noDists + 4] = m_f[cellId * m_noDists + 6];
////        m_fold[cellId * m_noDists + 7] = m_f[cellId * m_noDists + 5];
////      }
////      // +y
////      if(grid().neighbor(cellId, 3) == INVALID_CELLID) {
////        //-y = +y
////        m_fold[cellId * m_noDists + 2] = m_f[cellId * m_noDists + 3];
////        // +x-y = -x+y
////        m_fold[cellId * m_noDists + 5] = m_f[cellId * m_noDists + 7] + 1.0 / 6.0 * (0.1); // todo:testing
////        // +x-y = -x+y
////        m_fold[cellId * m_noDists + 6] = m_f[cellId * m_noDists + 4] - 1.0 / 6.0 * (0.1); // todo:testing
////      }
//    }
//  }
}

template <Debug_Level DEBUG_LEVEL>
void LBMSolver<DEBUG_LEVEL>::transferGrid(const GridInterface& grid) {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMInit]);
  m_dim = grid.dim();

  cerr0 << "Transferring " + std::to_string(m_dim) + "D Grid to LBM solver" << std::endl;
  logger << "Transferring " + std::to_string(m_dim) + "D Grid to LBM solver" << std::endl;


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
}

template <Debug_Level DEBUG_LEVEL>
constexpr auto LBMSolver<DEBUG_LEVEL>::isThermal() -> GBool {
  return false;
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
auto LBMSolver<DEBUG_LEVEL>::convergence(const GInt var) const -> GDouble {
  GDouble conv = 0.0;
  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    conv += std::abs(m_vars[cellId * m_noVars + var] - m_varsold[cellId * m_noVars + var]);
  }
  return conv;
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void LBMSolver<DEBUG_LEVEL>::setupPeriodicConnection() {
  std::unordered_multimap<GDouble, GInt> coordToCellIds;
  std::set<GDouble>                      keys;
  std::set<GInt>                      cellList;

  for(GInt cellId = 0; cellId < noCells(); ++cellId) {
    if(grid().property(cellId, CellProperties::bndry) && grid().property(cellId, CellProperties::leaf)) {
      if(grid().neighbor(cellId, 0) == INVALID_CELLID || grid().neighbor(cellId, 1) == INVALID_CELLID) {
        coordToCellIds.emplace(grid<NDIM>()->center(cellId, 1), cellId);
        keys.emplace(grid<NDIM>()->center(cellId, 1));
        cellList.emplace(cellId);
      }
    }
  }
  for(const auto& p : keys) {
    if(coordToCellIds.count(p) == 2) {
      auto search = coordToCellIds.equal_range(p);
      for(auto it = search.first; it != search.second; ++it) {
        if(grid().neighbor(it->second, 0) == INVALID_CELLID) {
          const GInt neighborA                 = it->second;
          const GInt neighborB                 = (++it)->second;
          grid<NDIM>()->neighbor(neighborA, 0) = neighborB;
          grid<NDIM>()->neighbor(neighborB, 1) = neighborA;
        } else {
          const GInt neighborA                 = it->second;
          const GInt neighborB                 = (++it)->second;
          grid<NDIM>()->neighbor(neighborA, 1) = neighborB;
          grid<NDIM>()->neighbor(neighborB, 0) = neighborA;
        }

      }
    } else {
      TERMM(-1, "Invalid result!" + std::to_string(p));
    }
  }

  //todo: move this periodic code before diagonal detection
  for(const auto cellId:cellList){
    const GInt nghbrmX = grid<NDIM>()->neighbor(cellId, 0);
    const GInt nghbrpX = grid<NDIM>()->neighbor(cellId, 1);

    // +x+y
    const GInt nghbrpXpY                                 = nghbrpX != INVALID_CELLID ? grid<NDIM>()->neighbor(nghbrpX, 3) : -1;
    grid<NDIM>()->neighbor(cellId, cartesian::maxNoNghbrs<NDIM>()) = nghbrpXpY;

    // +x-y
    const GInt nghbrpXmY                                 = nghbrpX != INVALID_CELLID ? grid<NDIM>()->neighbor(nghbrpX, 2) : -1;
    grid<NDIM>()->neighbor(cellId, cartesian::maxNoNghbrs<NDIM>() + 1) = nghbrpXmY;

    // -x-y
    const GInt nghbrmXmY                                 = nghbrmX != INVALID_CELLID ? grid<NDIM>()->neighbor(nghbrmX, 2) : -1;
    grid<NDIM>()->neighbor(cellId, cartesian::maxNoNghbrs<NDIM>() + 2) = nghbrmXmY;

    // -x+y
    const GInt nghbrmXpY                                 = nghbrmX != INVALID_CELLID ? grid<NDIM>()->neighbor(nghbrmX, 3) : -1;
    grid<NDIM>()->neighbor(cellId, cartesian::maxNoNghbrs<NDIM>() + 3) = nghbrmXpY;
  }
}

template class LBMSolver<Debug_Level::no_debug>;
template class LBMSolver<Debug_Level::min_debug>;
template class LBMSolver<Debug_Level::debug>;
template class LBMSolver<Debug_Level::more_debug>;
template class LBMSolver<Debug_Level::max_debug>;