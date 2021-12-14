#include "lbm_solver.h"
#include "analytical_solutions.h"
#include "lbm_equilibrium_func.h"

#include <set>

using namespace std;

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::init(int argc, GChar** argv) {
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

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::init(int argc, GChar** argv, GString config_file) {
  setConfiguration(config_file);
  init(argc, argv);
  logger << NDIM << "D LBM Solver started ||>" << endl;
  cout << NDIM << "D LBM Solver started ||>" << endl;
  Configuration::load("solver");
  POST::setConfAccessor(Configuration::getAccessor("postprocessing"));
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMInit]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::initBenchmark(int argc, GChar** argv) {
  m_benchmark = true;
  init(argc, argv);

  logger << "Setting up benchmarking!" << endl;

  TERMM(-1, "Not implemented!");
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::initTimers() {
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMSolverTotal], "Total run time of the LBM Solver.", TimeKeeper[Timers::timertotal]);
  RECORD_TIMER_START(TimeKeeper[Timers::LBMSolverTotal]);

  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMInit], "Initialization of the LBM solver!", TimeKeeper[Timers::LBMSolverTotal]);
  RECORD_TIMER_START(TimeKeeper[Timers::LBMInit]);

  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMMainLoop], "Main Loop of the LBM solver!", TimeKeeper[Timers::LBMSolverTotal]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMCalc], "Computation", TimeKeeper[Timers::LBMMainLoop]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMColl], "Collision", TimeKeeper[Timers::LBMCalc]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMProp], "Propagation", TimeKeeper[Timers::LBMCalc]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMBnd], "Boundary condition", TimeKeeper[Timers::LBMCalc]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMEq], "Equilibrium distribution", TimeKeeper[Timers::LBMCalc]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMMacro], "Computing macroscopic values", TimeKeeper[Timers::LBMCalc]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMForce], "Computing forcing", TimeKeeper[Timers::LBMCalc]);


  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMPost], "Postprocessing", TimeKeeper[Timers::LBMMainLoop]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::LBMIo], "IO", TimeKeeper[Timers::LBMMainLoop]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::loadConfiguration() {
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
  if(EQ == LBEquation::Navier_Stokes || EQ == LBEquation::Navier_Stokes_Poisson) {
    m_re = required_config_value<GDouble>("reynoldsnumber");
  } else {
    m_re = 0;
  }


  /// todo: 1.73205080756887729352??? (F1BCS)
  //  m_nu    = m_ma / 1.73205080756887729352 / m_re * m_referenceLength; //* (FFPOW2[maxLevel() - a_level(pCellId)]);
  //  m_omega = 2.0 / (1.0 + 6.0 * m_nu);

  m_omega = 1.0 / m_relaxTime;
  m_nu    = (2.0 * m_relaxTime - 1) / 6.0;
  m_refU  = m_re * m_nu / m_refLength;
  m_dt    = 1.0 / (std::pow(size(), 1.0 / NDIM) - 1); // todo:test

  m_bndManager = std::make_unique<LBMBndManager<DEBUG_LEVEL, LBTYPE, EQ>>();

  // todo: doesn't work??
  //   std::function<const Surface<NDIM>&(GInt)> _bndrySurface = [this](const GInt id) { return bndrySurface(id); };
  using namespace placeholders;
  std::function<const Surface<DEBUG_LEVEL, NDIM>&(GString)> _bndrySurface = std::bind(&LBMSolver::bndrySurface, this, _1);
  m_bndManager->setupBndryCnds(required_config_value<json>("boundary"), _bndrySurface);


  cerr0 << "<<<<<<<<<<<<>>>>>>>>>>>>>" << std::endl;
  cerr0 << "LBM Type "
        << "BGK" << std::endl; // todo: fix me
  cerr0 << "LBM Model " << METH::m_name << std::endl;
  cerr0 << "Equation " << LBEquationName[static_cast<GInt>(EQ)] << std::endl;
  cerr0 << "No. of variables " << std::to_string(NVARS) << std::endl;
  cerr0 << "No. Leaf cells: " << noLeafCells() << std::endl;
  cerr0 << "No. Bnd cells: " << noBndCells() << std::endl;
  cerr0 << "Relaxation Time: " << m_relaxTime << std::endl;
  cerr0 << "Timestep: " << m_dt << std::endl;
  cerr0 << "Omega: " << m_omega << std::endl;
  cerr0 << "Reynolds Number: " << m_re << std::endl;
  cerr0 << "Viscosity: " << m_nu << std::endl;
  cerr0 << "Reference V: " << m_refU << std::endl;
  cerr0 << "+++++++++++++++++++++++++" << std::endl;
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::allocateMemory() {
  // allocate memory
  m_f.resize(grid().totalSize() * NDIST);
  m_feq.resize(grid().totalSize() * NDIST);
  m_fold.resize(grid().totalSize() * NDIST);
  m_vars.resize(grid().totalSize() * NVARS);
  m_varsold.resize(grid().totalSize() * NVARS);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
auto LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::run() -> GInt {
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

    // also write an output if it's the last time step or if the solution has converged
    output(m_timeStep == noTimesteps - 1 || converged);
  }

  if(has_config_value("analyticalSolution")) {
    compareToAnalyticalResult();
  }

  executePostprocess(pp::HOOK::ATEND);

  logger << "LBM Solver finished <||" << endl;
  cout << "LBM Solver finished <||" << endl;
  unusedConfigValues();
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMMainLoop]);
  return 0;
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
auto LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::convergenceCondition() -> GBool {
  if(m_timeStep > 0 && m_timeStep % m_outputInfoInterval == 0) {
    cerr0 << m_timeStep << ": ";
    logger << m_timeStep << ": ";

    std::array<GDouble, NVARS> conv;
    for(GInt vid = 0; vid < NVARS; ++vid) {
      conv[vid] = sumAbsDiff(vid);
      cerr0 << "d" << VAR::varStr(vid) << "=" << conv[vid] << " ";
      logger << "d" << VAR::varStr(vid) << "=" << conv[vid] << " ";
    }
    cerr0 << std::endl;
    logger << std::endl;

    const GDouble maxConv         = *std::max_element(conv.begin(), conv.end());
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

// todo: refactor
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::initialCondition() {
  // init to zero:
  for(GInt cellId = 0; cellId < allCells(); ++cellId) {
    for(const auto dir : VAR::velocities()) {
      m_vars[cellId * NVARS + dir]    = 0;
      m_varsold[cellId * NVARS + dir] = 0;
    }
    if(EQ == LBEquation::Navier_Stokes) {
      rho(cellId) = 1.0;
    } else if(EQ == LBEquation::Poisson || EQ == LBEquation::Navier_Stokes_Poisson) {
      electricPotential(cellId) = 0.0;
    }
  }
  initBndryValues();

  for(GInt cellId = 0; cellId < allCells(); ++cellId) {
    if(EQ == LBEquation::Navier_Stokes) {
      // assuming initial zero velocity and density 1
      for(GInt dist = 0; dist < NDIST; ++dist) {
        feq(cellId, dist)  = METH::m_weights[dist];
        f(cellId, dist)    = feq(cellId, dist);
        fold(cellId, dist) = feq(cellId, dist);
      }
    }
    if(EQ == LBEquation::Poisson) {
      for(GInt dist = 0; dist < NDIST - 1; ++dist) {
        feq(cellId, dist)  = METH::m_weights[dist] * electricPotential(cellId);
        f(cellId, dist)    = feq(cellId, dist);
        fold(cellId, dist) = feq(cellId, dist);
      }
      const GDouble centerFeq = (METH::m_weights[NDIST - 1] - 1.0) * electricPotential(cellId);
      feq(cellId, NDIST - 1)  = centerFeq;
      f(cellId, NDIST - 1)    = centerFeq;
      fold(cellId, NDIST - 1) = centerFeq;
    }
  }
}

/// Set the initial values on the boundaries
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::initBndryValues() {
  using namespace std::placeholders;
  std::function<GDouble&(GInt, GInt)> _v = std::bind(&LBMSolver::vars, this, _1, _2);

  m_bndManager->initCndBnd(_v);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::timeStep() {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMCalc]);
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
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMCalc]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::output(const GBool forced, const GString& postfix) {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMIo]);

  // todo: make the output values settable
  if(m_diverged) {
    // output the values before updating the macroscopic values in case the solution diverged
    output(true, "bdiv");
  }
  //  updateMacroscopicValues();
  if((m_timeStep > 0 && m_timeStep % m_outputSolutionInterval == 0) || forced) {
    std::vector<IOIndex>              index;
    std::vector<std::vector<GString>> values;

    // only output leaf cells (i.e. cells without children)
    std::function<GBool(GInt)> isLeaf = [&](GInt cellId) { return noChildren(cellId) == 0; };

    if(EQ != LBEquation::Poisson) {
      std::array<std::vector<GDouble>, NDIM> tmpVel;
      std::vector<GDouble>                   tmpRho;


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
        index.emplace_back(IOIndex{static_cast<GString>(VAR::VELSTR[dir]), "float32"});
        values.emplace_back(toStringVector(tmpVel[dir], size()));
      }
      // todo: fix type
      index.emplace_back(IOIndex{"rho", "float32"});
      values.emplace_back(toStringVector(tmpRho, size()));
    }

    if(EQ == LBEquation::Poisson || EQ == LBEquation::Navier_Stokes_Poisson) {
      std::vector<GDouble> tmpV;
      tmpV.resize(size());

      for(GInt cellId = 0; cellId < size(); ++cellId) {
        tmpV[cellId] = electricPotential(cellId);
      }
      // todo: fix type
      index.emplace_back(IOIndex{"V", "float32"});
      values.emplace_back(toStringVector(tmpV, size()));
    }

    VTK::ASCII::writePoints<NDIM>("test_" + std::to_string(m_timeStep) + postfix, size(), center(), index, values, isLeaf);
  }
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMIo]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::compareToAnalyticalResult() {
  // todo: fix for 3D
  //  for an analytical testcase the value "analyticalSolution" needs to be defined
  const auto analyticalSolutionName = required_config_value<GString>("analyticalSolution");
  // from the analyticalSolutionName we can get a functional representation of the solution
  const auto anaSolution = analytical::getAnalyticalSolution<NDIM>(analyticalSolutionName);

  // determine the error
  GDouble maxError = 0;

  // exclude cells from error calculation, for example, the boundary surfaces
  std::vector<std::vector<GInt>> m_excludedCells;
  if(has_config_value("analyticalSolutionExcludeSurface")) {
    const auto excludedSurfaces = required_config_value<std::vector<GString>>("analyticalSolutionExcludeSurface");
    for(const auto& exclSurfName : excludedSurfaces) {
      m_excludedCells.emplace_back(bndrySurface(exclSurfName).getCellList());
    }
  }

  GDouble sumError      = 0;
  GDouble sumErrorSq    = 0;
  GDouble sumSolution   = 0;
  GDouble sumSolutionSq = 0;
  for(GInt cellId = 0; cellId < noInternalCells(); ++cellId) {
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
    VectorD<NDIM> v(&velocity(cellId, 0));

    const Point<NDIM>   shiftedCenter = shiftCenter(center(cellId));
    const VectorD<NDIM> solution      = anaSolution(shiftedCenter);
    const GDouble       delta         = (v - solution).norm();
    sumError += delta;
    sumErrorSq += gcem::pow(delta, 2);
    sumSolution += solution.norm();
    sumSolutionSq += gcem::pow(solution.norm(), 2);
    maxError = std::max(delta, maxError);
  }
  const GDouble L2error = (gcem::sqrt(sumErrorSq) / gcem::sqrt(sumSolutionSq)) / gcem::pow(noInternalCells(), 1.0 / NDIM);

  // if that compare with the maximum expected error to give a pass/no-pass result
  const GDouble maximumExpectedError    = opt_config_value("errorMax", 1.0);
  const GDouble maximumExpectedErrorL2  = opt_config_value("errorL2", 1.0);
  const GDouble maximumExpectedErrorGRE = opt_config_value("errorGRE", 1.0);

  const GDouble gre = sumError / sumSolution;

  cerr0 << "Comparing to analytical result " << analyticalSolutionName << std::endl;
  logger << "Comparing to analytical result " << analyticalSolutionName << std::endl;
  cerr0 << "max. Error: " << maxError << std::endl;
  logger << "max. Error: " << maxError << std::endl;
  cerr0 << "avg. L2: " << L2error << std::endl;
  logger << "avg. L2: " << L2error << std::endl;
  cerr0 << "global relative error: " << gre << std::endl;
  logger << "global relative error: " << L2error << std::endl;

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
  if(maximumExpectedErrorGRE < gre) {
    failedErrorCheck = true;
    cerr0 << "Error bounds for the global relative error have failed!!! ( < " << maximumExpectedErrorGRE << ")" << std::endl;
    logger << "Error bounds for the global relative error have failed!!! ( < " << maximumExpectedErrorGRE << ")" << std::endl;
  }

  if(failedErrorCheck) {
    TERMM(-1, "Analytical testcase failed");
  }
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
auto LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::shiftCenter(const Point<NDIM>& centerPoint) const -> Point<NDIM> {
  const auto analyticalSolutionName = required_config_value<GString>("analyticalSolution");

  Point<NDIM> shiftedCenter = centerPoint;
  if(analyticalSolutionName == analytical::getStr(analytical::ANALYTICAL_CASE_INDEX::poissonCHAI08_1)) {
    if(NDIM > 1) {
      TERMM(-1, "FIX ME");
    }

    // align points equidistantly between (0, 1)
    const GDouble actualExtent    = std::abs(center(0, 0) - center(size() - 1, 0));
    const GDouble correctedExtent = 1.0;
    const GInt    pos             = std::floor(shiftedCenter[0] / (actualExtent / (size() - 1)));
    shiftedCenter[0]              = pos * correctedExtent / (size() - 1);
  }

  return shiftedCenter;
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::currToOldVars() {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMMacro]);
  std::copy(m_vars.begin(), m_vars.end(), m_varsold.begin());
  //  std::copy(m_f.begin(), m_f.end(), m_fold.begin());
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMMacro]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::updateMacroscopicValues() {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMMacro]);

#ifdef _OPENMP
#pragma omp parallel default(none)
  {
#endif

#ifdef _OPENMP
#pragma omp for
#endif
    for(GInt cellId = 0; cellId < allCells(); ++cellId) {
      // todo:skip non-leaf cells!
      if(EQ == LBEquation::Navier_Stokes || EQ == LBEquation::Navier_Stokes_Poisson) {
        rho(cellId) = std::accumulate((m_fold.begin() + cellId * NDIST), (m_fold.begin() + (cellId + 1) * NDIST), 0.0);

        for(GInt dir = 0; dir < NDIM; ++dir) {
          velocity(cellId, dir) = 0;
          for(GInt dist = 0; dist < NDIST - 1; ++dist) {
            velocity(cellId, dir) += METH::m_dirs[dist][dir] * fold(cellId, dist);
          }
          velocity(cellId, dir) /= rho(cellId);
        }
      }

      if(EQ == LBEquation::Poisson || EQ == LBEquation::Navier_Stokes_Poisson) {
        electricPotential(cellId) = 1.0 / (1.0 - METH::m_weights[NDIST - 1])
                                    * std::accumulate((m_fold.begin() + cellId * NDIST), (m_fold.begin() + (cellId + 1) * NDIST - 1), 0.0);
      }
    }
#ifdef _OPENMP
  }
#endif

  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMMacro]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::calcEquilibriumMoments() {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMEq]);

#ifdef _OPENMP
#pragma omp parallel for default(none)
#endif
  for(GInt cellId = 0; cellId < allCells(); ++cellId) {
    if(EQ != LBEquation::Poisson) {
      std::array<GDouble, NDIM> velos; // NOLINT(cppcoreguidelines-pro-type-member-init)
      for(GInt dir = 0; dir < NDIM; ++dir) {
        velos[dir] = velocity(cellId, dir);
      }

      //    const GDouble squaredV = gcem::pow(velocity<NDIM>(cellId, 0), 2) + gcem::pow(velocity<NDIM>(cellId, 1), 2);
      const GDouble density = rho(cellId);

      for(GInt dist = 0; dist < NDIST; ++dist) {
        feq(cellId, dist) = 1.0;
        for(GInt dir = 0; dir < NDIM; ++dir) {
          feq(cellId, dist) += 3 * METH::m_dirs[dist][dir] * velos[dir];
        }
        feq(cellId, dist) *= METH::m_weights[dist] * density;
        //      feq(cellId, dist) = METH::m_weights[dist] * density * (1 + 3 * (u * METH::m_dirs[dist][0] + v * METH::m_dirs[dist][1]));
        //      const GDouble cu             = velocity<NDIM>(cellId, 0) * cx[dist] + velocity<NDIM>(cellId, 1) * cy[dist];
        //      m_feq[cellId * NDIST + dist] = LBMethod<LBTYPE>::m_weights[dist] * rho<NDIM>(cellId) * (1 + 3 * cu + 4.5 * cu * cu - 1.5 *
        //      squaredV);
      }
    }
    if(EQ == LBEquation::Poisson) {
      eq::poisson<LBTYPE>(feq(cellId), electricPotential(cellId));
    }
  }
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMEq]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::collisionStep() {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMColl]);

  if(DEBUG_LEVEL > Debug_Level::debug) {
    if(hasNAN(m_fold) >= 0) {
      TERMM(-1, "NAN within fold");
    }
    if(hasNAN(m_feq) >= 0) {
      TERMM(-1, "NAN within feq");
    }
    if(hasNAN(m_vars) >= 0) {
      TERMM(-1, "NAN within vars");
    }
  }

  for(GInt cellId = 0; cellId < allCells(); ++cellId) {
    for(GInt dist = 0; dist < NDIST; ++dist) {
      f(cellId, dist) = (1 - m_omega) * fold(cellId, dist) + m_omega * feq(cellId, dist);
      if(EQ == LBEquation::Poisson && dist != NDIST - 1) {
        // todo: make settable
        static constexpr GDouble k           = 27.79;
        static const GDouble     diffusivity = METH::m_poissonAlpha * gcem::pow(m_speedOfSound, 2) * (0.5 - m_relaxTime) * m_dt;
        f(cellId, dist) += m_dt * diffusivity * METH::m_poissonWeights[dist] * k * k * electricPotential(cellId);
      }
    }
  }

  if(DEBUG_LEVEL > Debug_Level::debug) {
    if(hasNAN(m_f) >= 0) {
      TERMM(-1, "NAN within f");
    }
  }

  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMColl]);
}

// todo:refactor
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::forcing() {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMForce]);

  static GBool info        = true;
  const auto   forcingType = opt_config_value<GString>("forcing", "");

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
    const GDouble inletPressure = 1.0115985357; // 0.01 (D2Q9)
                                                //    const GDouble inletPressure = 1.00055; // 0.01


    // set Outlet forcing
    for(const GInt inletCellId : inlet.getCellList()) {
      // cell to the inside of the boundary used for extrapolation
      const GInt           valueCellId = grid().neighbor(inletCellId, 1);
      const VectorD<NDIM>& centerInlet = grid().center(valueCellId);
      for(const GInt outletCellId : outlet.getCellList()) {
        const VectorD<NDIM>& centerOutlet = grid().center(outletCellId);
        // todo: improve this since this searches for the outlet cellid every time step!
        if(std::abs(centerInlet[1] - centerOutlet[1]) < GDoubleEps) {
          for(GInt dist = 0; dist < LBMethod<LBTYPE>::m_noDists; ++dist) {
            // w(k)*(rho_inlet+ 3*(cx(k)*u(NX-1,:)+cy(k)*v(NX-1,:)))+(f(NX-1,:,k)-feq(NX-1,:,k))
            // recalculate f on the outlet cell using a given outlet pressure
            //            f(outletCellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * outletPressure
            //                                        * (1.0
            //                                           + 3.0 * (vars(valueCellId, VAR::velocity(0)) * LBMethod<LBTYPE>::m_dirs[dist][0]
            //                                                    + vars(valueCellId, VAR::velocity(1)) *
            //                                                    LBMethod<LBTYPE>::m_dirs[dist][1]))
            //                                    + f(valueCellId, dist) - feq(valueCellId, dist);
            f(outletCellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * outletPressure
                                        * (1
                                           + 3
                                                 * (vars(valueCellId, VAR::velocity(0)) * LBMethod<LBTYPE>::m_dirs[dist][0]
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
                                                * (vars(valueCellId, VAR::velocity(0)) * LBMethod<LBTYPE>::m_dirs[dist][0]
                                                   + 0 * LBMethod<LBTYPE>::m_dirs[dist][1]))
                                   + f(valueCellId, dist) - feq(valueCellId, dist);
          }
        }
      }
    }
  }
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMForce]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::prePropBoundaryCnd() {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMBnd]);

  // todo: replace by lambda function
  using namespace std::placeholders;
  std::function<GDouble&(GInt, GInt)> _fold = std::bind(&LBMSolver::fold, this, _1, _2);
  std::function<GDouble&(GInt, GInt)> _f    = std::bind(&LBMSolver::f, this, _1, _2);
  std::function<GDouble&(GInt, GInt)> _v    = std::bind(&LBMSolver::vars, this, _1, _2);

  m_bndManager->preApply(_f, _fold, _v);
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMBnd]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::propagationStep() {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMProp]);

  if(DEBUG_LEVEL > Debug_Level::debug) {
    if(hasNAN(m_fold) >= 0) {
      TERMM(-1, "NAN within fold");
    }
    if(hasNAN(m_f) >= 0) {
      TERMM(-1, "NAN within f");
    }
  }

#ifdef _OPENMP
#pragma omp parallel for default(none)
#endif
  for(GInt cellId = 0; cellId < noInternalCells(); ++cellId) {
    for(GInt dist = 0; dist < NDIST - 1; ++dist) {
      const GInt neighborId = m_grid->neighbor(cellId, dist);
      if(neighborId != INVALID_CELLID) {
        fold(neighborId, dist) = f(cellId, dist);
      }
    }
    fold(cellId, NDIST - 1) = f(cellId, NDIST - 1);
  }
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMProp]);
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::boundaryCnd() {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMBnd]);

  // todo: replace by lambda function
  using namespace std::placeholders;
  std::function<GDouble&(GInt, GInt)> _fold = std::bind(&LBMSolver::fold, this, _1, _2);
  std::function<GDouble&(GInt, GInt)> _f    = std::bind(&LBMSolver::f, this, _1, _2);
  std::function<GDouble&(GInt, GInt)> _v    = std::bind(&LBMSolver::vars, this, _1, _2);

  m_bndManager->apply(_f, _fold, _v);
  RECORD_TIMER_STOP(TimeKeeper[Timers::LBMBnd]);
}

// todo: simplify (move somewhere else?!)
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
void LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::transferGrid(const GridInterface& grid) {
  RECORD_TIMER_START(TimeKeeper[Timers::LBMInit]);
  cerr0 << "Transferring " + std::to_string(NDIM) + "D Grid to LBM solver" << std::endl;
  logger << "Transferring " + std::to_string(NDIM) + "D Grid to LBM solver" << std::endl;

  if(grid.dim() != NDIM) {
    TERMM(-1, "Invalid configuration the grid dimensionality is not matching!");
  }


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

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
constexpr auto LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::isThermal() -> GBool {
  return false;
}

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
auto LBMSolver<DEBUG_LEVEL, LBTYPE, EQ>::sumAbsDiff(const GInt var) const -> GDouble {
  GDouble conv = 0.0;
  for(GInt cellId = 0; cellId < noInternalCells(); ++cellId) {
    conv += std::abs(m_vars[cellId * NVARS + var] - m_varsold[cellId * NVARS + var]);
  }
  return conv;
}

// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q15>;
// template class LBMSolver<Debug_Level::min_debug, LBMethodType::D3Q15>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q15>;
// template class LBMSolver<Debug_Level::more_debug, LBMethodType::D3Q15>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q15>;
