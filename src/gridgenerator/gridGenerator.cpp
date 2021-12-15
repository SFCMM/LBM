#include <iostream>
#include <mpi.h>
#include <utility>

#include <sfcmm_common.h>

#include "cartesiangrid_generation.h"
#include "config.h"
#include "geometry.h"
#include "gridGenerator.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

template <Debug_Level DEBUG_LEVEL>
void GridGenerator<DEBUG_LEVEL>::init(int argc, GChar** argv) {
  m_exe = argv[0]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

#ifndef GRIDGEN_SINGLE_FILE_LOG
  logger.open("gridgen_log" + std::to_string(m_domainId), false, argc, argv, MPI_COMM_WORLD);
#else
  if(DEBUG_LEVEL < Debug_Level::max_debug) {
    logger.open("gridgen_log", true, argc, argv, MPI_COMM_WORLD);
  } else {
    logger.open("gridgen_log", false, argc, argv, MPI_COMM_WORLD);
  }
#endif
  logger.setMinFlushSize(LOG_MIN_FLUSH_SIZE);

  initTimers();
}

template <Debug_Level DEBUG_LEVEL>
void GridGenerator<DEBUG_LEVEL>::init(int argc, GChar** argv, GString config_file) {
  setConfiguration(config_file);
  init(argc, argv);
}

template <Debug_Level DEBUG_LEVEL>
void GridGenerator<DEBUG_LEVEL>::initBenchmark(int argc, GChar** argv) {
  m_benchmark = true;

  m_dim              = 3;
  m_maxNoCells       = 100000;
  m_partitionLvl     = 3;
  m_uniformLvl       = 5;
  m_maxRefinementLvl = m_uniformLvl;

  init(argc, argv);
}

template <Debug_Level DEBUG_LEVEL>
void GridGenerator<DEBUG_LEVEL>::initTimers() {
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridGeneratorTotal], "Total run time of the grid generator", TimeKeeper[Timers::timertotal]);
  RECORD_TIMER_START(TimeKeeper[Timers::GridGeneratorTotal]);

  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::Init], "Init", TimeKeeper[Timers::GridGeneratorTotal]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridGeneration], "Create the grid.", TimeKeeper[Timers::GridGeneratorTotal]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridInit], "Init grid.", TimeKeeper[Timers::GridGeneration]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridPart], "Partitioning grid generation.", TimeKeeper[Timers::GridGeneration]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridUniform], "Uniform grid generation.", TimeKeeper[Timers::GridGeneration]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridRefinement], "Grid refinement.", TimeKeeper[Timers::GridGeneration]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridIo], "Grid IO.", TimeKeeper[Timers::GridGeneratorTotal]);
  NEW_TIMER_NOCREATE(TimeKeeper[Timers::IO], "IO", TimeKeeper[Timers::timertotal]);
}

template <Debug_Level DEBUG_LEVEL>
auto GridGenerator<DEBUG_LEVEL>::run() -> GInt {
  TimerProfiling runProfile("GridGenerator::run");
  PROFILE();
  RECORD_TIMER_START(TimeKeeper[Timers::Init]);
  logger << "Grid generator started ||>" << endl;
  cout << "Grid generator started ||>" << endl;
  loadConfiguration();
  RECORD_TIMER_STOP(TimeKeeper[Timers::Init]);

  switch(m_dim) {
    case 1:
      generateGrid<1>();
      break;
    case 2:
      generateGrid<2>();
      break;
    case 3:
      generateGrid<3>();
      break;
    case 4:
      generateGrid<4>();
      break;
    default:
      TERMM(-1, "Invalid number of dimensions 1-4.");
  }

  logger << "Grid generator finished <||" << endl;
  cout << "Grid generator finished <||" << endl;

  unusedConfigValues();

  return 0;
}

template <Debug_Level DEBUG_LEVEL>
void GridGenerator<DEBUG_LEVEL>::loadConfiguration() {
  RECORD_TIMER_START(TimeKeeper[Timers::IO]);
  RECORD_TIMER_START(TimeKeeper[Timers::GridIo]);

  if(!m_benchmark) {
    Configuration::load();
  } else {
    logger << "Setting up benchmarking!" << endl;
  }

  // 2. communicate configuration file to all other processes
  if(!MPI::isSerial()) {
    TERMM(-1, "Not implemented!");
  }

  // 3. load&check configuration values
  if(!m_benchmark) {
    m_dim        = required_config_value<GInt>("dim");
    m_maxNoCells = required_config_value<GInt>("maxNoCells");
  }

  // todo: unused
  m_dryRun = opt_config_value<GBool>("dry-run", m_dryRun);

  m_outputDir = getCWD() + "/" + opt_config_value<GString>("outputDir", m_outputDir);
  if(!isPath(m_outputDir, true)) {
    TERMM(-1, "Is not a valid output directory! " + m_outputDir);
  }

  // fix path
  if(m_outputDir.back() != '/') {
    m_outputDir += '/';
  }

  RECORD_TIMER_STOP(TimeKeeper[Timers::IO]);
  RECORD_TIMER_STOP(TimeKeeper[Timers::GridIo]);
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void GridGenerator<DEBUG_LEVEL>::generateGrid() {
  PROFILE();
  RECORD_TIMER_START(TimeKeeper[Timers::GridGeneration]);
  RECORD_TIMER_START(TimeKeeper[Timers::GridInit]);

  logger << "Generating a grid[" << NDIM << "D]" << endl;


  cout << SP1 << "Reading Grid definition" << endl;
  m_grid = std::make_unique<CartesianGridGen<DEBUG_LEVEL, NDIM>>(m_maxNoCells);
  if(!m_benchmark) {
    loadGridDefinition<NDIM>();
  } else {
    benchmarkSetup<NDIM>();
  }
  logger << SP2 << "+ maximum number of cells: " << m_maxNoCells << endl;
  cout << SP2 << "+ maximum number of cells: " << m_maxNoCells << endl;

  // todo: add function to define memory to be allocated
  // todo: add function to convert to appropriate memory size
  const GDouble memoryConsumptionKB = CartesianGridGen<DEBUG_LEVEL, NDIM>::memorySizePerCell() * m_maxNoCells / DKBIT;

  logger << SP2 << "+ local memory allocated: " << memoryConsumptionKB << "KB" << std::endl;
  cout << SP2 << "+ local memory allocated: " << memoryConsumptionKB << "KB" << std::endl;

  if(!MPI::isSerial()) {
    const GDouble globalMemory = memoryConsumptionKB * static_cast<GDouble>(MPI::globalNoDomains());
    logger << SP2 << "+ global memory allocated: " << globalMemory << "KB" << std::endl;
    cout << SP2 << "+ global memory allocated: " << globalMemory << "KB" << std::endl;
  }

  m_grid->setMaxLvl(m_maxRefinementLvl);
  // todo: allow setting the weighting method
  m_weightMethod = std::make_unique<WeightUniform>();


  logger << "\n";
  logger << SP2 << "+ m_center of gravity: " << strStreamify<NDIM>(m_grid->cog()).str() << "\n";
  logger << SP2 << "+ decisive direction: " << m_grid->largestDir() << "\n";
  logger << SP2 << "+ geometry extents: " << strStreamify<NDIM>(m_grid->lengthOfBoundingBox()).str() << "\n";
  logger << SP2 << "+ bounding box: " << m_grid->boundingBox().str() << endl;
  RECORD_TIMER_STOP(TimeKeeper[Timers::GridInit]);

  const std::function<GString()> strHighestLvl = [&]() { return to_string(m_grid->currentHighestLvl()); };
  logger.addAttribute({"level", strHighestLvl});

  // create partitioning grid first, which is done without MPI parallelization
  gridGen<NDIM>().createPartitioningGrid(m_partitionLvl);

  if(!MPI::isSerial()) {
    // todo:implement
    TERMM(-1, "Not implemented");
    // m_grid->setupMPIComm();
  }

  gridGen<NDIM>().uniformRefineGrid(m_uniformLvl);

  RECORD_TIMER_START(TimeKeeper[Timers::GridRefinement]);
  for(GInt refinedLvl = m_uniformLvl; refinedLvl < m_maxRefinementLvl; ++refinedLvl) {
    GInt noCellsToRefine = gridGen<NDIM>().markBndryCells();
    gridGen<NDIM>().refineMarkedCells(noCellsToRefine);
    logger.updateAttributes();
  }
  RECORD_TIMER_STOP(TimeKeeper[Timers::GridRefinement]);
  RECORD_TIMER_STOP(TimeKeeper[Timers::GridGeneration]);
  logger.eraseAttribute("level");

  RECORD_TIMER_START(TimeKeeper[Timers::IO]);
  RECORD_TIMER_START(TimeKeeper[Timers::GridIo]);
  m_grid->save(m_outputDir + m_outGridFilename, m_gridOutConfig);
  RECORD_TIMER_STOP(TimeKeeper[Timers::IO]);
  RECORD_TIMER_STOP(TimeKeeper[Timers::GridIo]);
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void GridGenerator<DEBUG_LEVEL>::loadGridDefinition() {
  RECORD_TIMER_START(TimeKeeper[Timers::IO]);
  RECORD_TIMER_START(TimeKeeper[Timers::GridIo]);
  if(m_benchmark) {
    return;
  }

  m_partitionLvl = required_config_value<GInt>("partitionLevel");
  m_uniformLvl   = required_config_value<GInt>("uniformLevel");
  if(m_partitionLvl > m_uniformLvl) {
    TERMM(-1, "Invalid definition of grid level partitionLevel >= uniformLevel");
  }

  m_maxRefinementLvl = m_uniformLvl;
  m_maxRefinementLvl = opt_config_value<GInt>("maxRfnmtLvl", m_maxRefinementLvl);
  if(m_maxRefinementLvl < m_uniformLvl) {
    TERMM(-1, "Invalid definition of grid level uniformLevel >= maxRfnmtLvl");
  }

  // todo: unused
  m_maxNoOffsprings = opt_config_value<GInt>("maxNoOffsprings", m_maxNoOffsprings);

  m_outGridFilename = opt_config_value<GString>("gridFileName", m_outGridFilename);

  json defaultGridOutConfig = {{"format", "ASCII"}, {"cellFilter", "highestLvl"}, {"type", "points"}};
  m_gridOutConfig           = opt_config_value<json>("gridOutput", defaultGridOutConfig);

  cout << SP1 << "Reading Geometry" << endl;
  m_geometry = std::make_shared<GeometryManager<DEBUG_LEVEL, NDIM>>(MPI_COMM_WORLD);

  json defaultGeometry = {"cube", {{"center", {0.0, 0.0, 0.0}}, {"length", 1}}};
  m_geometryConfig     = opt_config_value<json>("geometry", defaultGeometry);
  m_geometry->setup(m_geometryConfig);
  m_grid->setGeometryManager(m_geometry);

  if(m_geometry->noObjects() == 0 || has_config_value("boundingBox")) {
    auto               tmp = opt_config_value<std::vector<GDouble>>("boundingBox", DEFAULT_BOUNDINGBOX.at(m_dim - 1));
    BoundingBoxDynamic tmpBB;
    tmpBB.init(m_dim);
    for(GInt dir = 0; dir < m_dim; ++dir) {
      tmpBB.min(dir) = tmp[2 * dir];
      tmpBB.max(dir) = tmp[2 * dir + 1];
    }

    m_grid->setBoundingBox(*static_cast<BoundingBoxInterface*>(static_cast<void*>(&tmpBB)));
  } else {
    auto tmp = m_geometry->getBoundingBox();
    m_grid->setBoundingBox(*static_cast<BoundingBoxInterface*>(static_cast<void*>(&tmp)));
  }


  RECORD_TIMER_STOP(TimeKeeper[Timers::IO]);
  RECORD_TIMER_STOP(TimeKeeper[Timers::GridIo]);
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void GridGenerator<DEBUG_LEVEL>::benchmarkSetup() {
  logger << "Setting up benchmarking grid!" << std::endl;
  m_geometry = std::make_shared<GeometryManager<DEBUG_LEVEL, NDIM>>(MPI_COMM_WORLD);

  json defaultGeometry = {{"cube", {{"type", "cube"}, {"center", {0.0, 0.0, 0.0}}, {"length", 1}}}};
  m_geometry->setup(defaultGeometry);
  m_grid->setGeometryManager(m_geometry);
  m_grid->setBoundingBox(m_geometry->getBoundingBox());
}


template class GridGenerator<Debug_Level::no_debug>;
template class GridGenerator<Debug_Level::debug>;
template class GridGenerator<Debug_Level::max_debug>;