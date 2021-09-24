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

#ifdef _OPENMP
  int provided = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
  MPI_Init(&argc, &argv);
#endif

  MPI_Comm_rank(MPI_COMM_WORLD, &m_domainId);
  MPI_Comm_size(MPI_COMM_WORLD, &m_noDomains);
  MPI::g_mpiInformation.init(m_domainId, m_noDomains);

  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  // Open cerr0 on MPI root
  if(MPI::isRoot()) {
    cerr0.rdbuf(std::cerr.rdbuf());
  } else {
    cerr0.rdbuf(&nullBuffer);
  }

#ifndef GRIDGEN_SINGLE_FILE_LOG
  logger.open("gridgen_log" + std::to_string(m_domainId), false, argc, argv, MPI_COMM_WORLD);
#else
  if(DEBUG_LEVEL < Debug_Level::more_debug) {
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
  m_configurationFileName = std::move(config_file);
  init(argc, argv);
}

template <Debug_Level DEBUG_LEVEL>
void GridGenerator<DEBUG_LEVEL>::initBenchmark(int argc, GChar** argv) {
  m_benchmark             = true;
  m_configurationFileName = "";

  m_dim              = 3;
  m_maxNoCells       = 100000;
  m_partitionLvl     = 3;
  m_uniformLvl       = 5;
  m_maxRefinementLvl = m_uniformLvl;

  init(argc, argv);
}

template <Debug_Level DEBUG_LEVEL>
void GridGenerator<DEBUG_LEVEL>::initTimers() {
  RESET_TIMERS();

  NEW_TIMER_GROUP_NOCREATE(TimeKeeper[Timers::AppGroup], "Application");
  NEW_TIMER_NOCREATE(TimeKeeper[Timers::timertotal], "Total", TimeKeeper[Timers::AppGroup]);
  RECORD_TIMER_START(TimeKeeper[Timers::timertotal]);

  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::Init], "Init", TimeKeeper[Timers::timertotal]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridGeneration], "Create the grid.", TimeKeeper[Timers::timertotal]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridInit], "Init grid.", TimeKeeper[Timers::GridGeneration]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridPart], "Partitioning grid generation.", TimeKeeper[Timers::GridGeneration]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridUniform], "Uniform grid generation.", TimeKeeper[Timers::GridGeneration]);
  NEW_SUB_TIMER_NOCREATE(TimeKeeper[Timers::GridRefinement], "Grid refinement.", TimeKeeper[Timers::GridGeneration]);
  NEW_TIMER_NOCREATE(TimeKeeper[Timers::IO], "IO", TimeKeeper[Timers::timertotal]);
}

template <Debug_Level DEBUG_LEVEL>
auto GridGenerator<DEBUG_LEVEL>::run() -> int {
  TimerProfiling runProfile("GridGenerator::run");
  PROFILE();
  RECORD_TIMER_START(TimeKeeper[Timers::Init]);
  startupInfo();
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

  RECORD_TIMER_STOP(TimeKeeper[Timers::timertotal]);
  STOP_ALL_RECORD_TIMERS();
  DISPLAY_ALL_TIMERS();

  MPI_Finalize();

  return 0;
}
template <Debug_Level DEBUG_LEVEL>
void GridGenerator<DEBUG_LEVEL>::startupInfo() {
  if(MPI::isRoot()) {
    cout << R"(    __  _______  __  _________     _     __)" << endl;
    cout << R"(   /  |/  / __ \/  |/  / ____/____(_)___/ /)" << endl;
    cout << R"(  / /|_/ / / / / /|_/ / / __/ ___/ / __  / )" << endl;
    cout << R"( / /  / / /_/ / /  / / /_/ / /  / / /_/ /  )" << endl;
    cout << R"(/_/  /_/\____/_/  /_/\____/_/  /_/\__,_/   )" << endl;
    cout << R"(                                           )" << endl;
    cout << "Start time:            " << dateString() << "\n"
         << "Number of ranks:       " << MPI::globalNoDomains() << "\n"
#ifdef _OPENMP
         << "Number of OMP threads: " << omp_get_max_threads() << "\n"
#endif
         << "Host (of rank 0):      " << hostString() << "\n"
         << "Working directory:     " << getCWD() << "\n"
         << "Executable:            " << m_exe << "\n"
         << endl;
  }
}

template <Debug_Level DEBUG_LEVEL>
void GridGenerator<DEBUG_LEVEL>::loadConfiguration() {
  RECORD_TIMER_START(TimeKeeper[Timers::IO]);

  if(!m_benchmark) {
    logger << "Loading configuration file [" << m_configurationFileName << "]" << endl;
  } else {
    logger << "Setting up benchmarking!" << endl;
  }

  // 1. open configuration file on root process
  if(MPI::isRoot() && isFile(m_configurationFileName)) {
    std::ifstream configFileStream(m_configurationFileName);
    configFileStream >> m_config;
    // put all available keys in map to keep track of usage
    for(const auto& element : m_config.items()) {
      m_configKeys.emplace(element.key(), false);
    }
    configFileStream.close();
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
  m_grid->createPartitioningGrid(m_partitionLvl);

  if(!MPI::isSerial()) {
    // todo:implement
    TERMM(-1, "Not implemented");
    // m_grid->setupMPIComm();
  }

  m_grid->uniformRefineGrid(m_uniformLvl);

  RECORD_TIMER_START(TimeKeeper[Timers::GridRefinement]);
  for(GInt refinedLvl = m_uniformLvl; refinedLvl < m_maxRefinementLvl; ++refinedLvl) {
    GInt noCellsToRefine = m_grid->markBndryCells();
    m_grid->refineMarkedCells(noCellsToRefine);
    logger.updateAttributes();
  }
  RECORD_TIMER_STOP(TimeKeeper[Timers::GridRefinement]);
  RECORD_TIMER_STOP(TimeKeeper[Timers::GridGeneration]);
  logger.eraseAttribute("level");

  RECORD_TIMER_START(TimeKeeper[Timers::IO]);
  m_grid->save(m_outputDir + m_outGridFilename, m_gridOutConfig);
  RECORD_TIMER_STOP(TimeKeeper[Timers::IO]);
}

template <Debug_Level DEBUG_LEVEL>
template <GInt NDIM>
void GridGenerator<DEBUG_LEVEL>::loadGridDefinition() {
  RECORD_TIMER_START(TimeKeeper[Timers::IO]);
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
      tmpBB.min[dir] = tmp[2 * dir];
      tmpBB.max[dir] = tmp[2 * dir + 1];
    }

    m_grid->setBoundingBox(*static_cast<BoundingBoxInterface*>(static_cast<void*>(&tmpBB)));
  } else {
    auto tmp = m_geometry->getBoundingBox();
    m_grid->setBoundingBox(*static_cast<BoundingBoxInterface*>(static_cast<void*>(&tmp)));
  }


  RECORD_TIMER_STOP(TimeKeeper[Timers::IO]);
}

template <Debug_Level DEBUG_LEVEL>
template <typename T>
auto GridGenerator<DEBUG_LEVEL>::required_config_value(const GString& key) -> T {
  // todo: check for types
  if(m_config.template contains(key)) {
    m_configKeys[key] = true;
    return static_cast<T>(m_config[key]);
  }
  TERMM(-1, "The required configuration value is missing: " + key);
}

template <Debug_Level DEBUG_LEVEL>
template <typename T>
auto GridGenerator<DEBUG_LEVEL>::opt_config_value(const GString& key, const T& defaultValue) -> T {
  // todo: check for types
  if(m_config.template contains(key)) {
    m_configKeys[key] = true;
    return static_cast<T>(m_config[key]);
  }
  return defaultValue;
}

template <Debug_Level DEBUG_LEVEL>
auto GridGenerator<DEBUG_LEVEL>::has_config_value(const GString& key) -> GBool {
  return m_config.template contains(key);
}


template <Debug_Level DEBUG_LEVEL>
void GridGenerator<DEBUG_LEVEL>::unusedConfigValues() {
  GInt i = 0;
  logger << "The following values in the configuration file are unused:" << endl;
  for(const auto& configKey : m_configKeys) {
    if(!configKey.second) {
      logger << "[" << ++i << "] " << configKey.first << "\n";
    }
  }
  logger << endl;
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
template class GridGenerator<Debug_Level::min_debug>;
template class GridGenerator<Debug_Level::debug>;
template class GridGenerator<Debug_Level::more_debug>;
template class GridGenerator<Debug_Level::max_debug>;