#include <cxxopts.hpp>
#include <mpi.h>

#include <sfcmm_common.h>
#include "config.h"
#include "interface/solver_interface.h"
#ifdef SOLVER_AVAILABLE
#include "gridgenerator/gridGenerator.h"
#include "lbm/lbm_solverExe.h"
#else
#include "gridGenerator.h"
#endif


#ifdef _OPENMP
#include <omp.h>
#endif

namespace internal_ {
/// Helper object to store some general application configuration options
template <SolverType Solver = SolverType::NONE>
class AppConfiguration {
 public:
  AppConfiguration(const int argc, GChar** argv, const GInt32 domainId, const GInt32 noDomains)
    : m_argc(argc), m_argv(argv), m_domainId(domainId), m_noDomains(noDomains) {}

  void init(const GInt debug) {
    switch(debug) {
      case static_cast<GInt>(Debug_Level::no_debug):
        init<Debug_Level::no_debug>();
        return;
      case static_cast<GInt>(Debug_Level::min_debug):
        init<Debug_Level::min_debug>();
        return;
      case static_cast<GInt>(Debug_Level::debug):
        init<Debug_Level::debug>();
        return;
      case static_cast<GInt>(Debug_Level::more_debug):
        init<Debug_Level::more_debug>();
        return;
      case static_cast<GInt>(Debug_Level::max_debug):
        [[fallthrough]];
      default:
        init<Debug_Level::max_debug>();
    }
  }

  auto run(const GInt debug) -> int {
    init(debug);
    return static_cast<int>(m_app->run());
  }

  void setConfigurationFile(const GString& configFile) { m_configurationFile = configFile; }
  void setBenchmark() { m_benchmark = true; }

  [[nodiscard]] auto toRun(const SolverType solver) const -> GBool {
    // open configuration file
    if(MPI::isRoot() && isFile(m_configurationFile)) {
      std::ifstream configFileStream(m_configurationFile);
      json          config;
      configFileStream >> config;
      // check for solvers for solver-> type
      const GString solverType = config["solver"]["type"];
      if(SOLVER_NAME[static_cast<GInt>(solver)] == solverType || SOLVER_NAMELC[static_cast<GInt>(solver)] == solverType) {
        return true;
      }
    }
    // todo: communicate the decision
    return false;
  }

  [[nodiscard]] auto grid() const -> const GridInterface& { return m_app->grid(); }

  void transferGrid(const GridInterface& grid, const GInt debug) {
    init(debug);
    m_app->transferGrid(grid);
  }

  void releaseMemory() { m_app.reset(nullptr); }

 private:
  template <Debug_Level DEBUG>
  void init() {
    if(m_init) {
      return;
    }

    switch(Solver) {
#ifdef SOLVER_AVAILABLE
      case SolverType::LBM:
        m_app = std::make_unique<LBMSolverExecutor<DEBUG>>(m_domainId, m_noDomains);
#endif
        break;
      case SolverType::NONE:
        [[fallthrough]];
      default:
        m_app = std::make_unique<GridGenerator<DEBUG>>(m_domainId, m_noDomains);
    }

    if(!m_benchmark) {
      m_app->init(m_argc, m_argv, m_configurationFile);
    } else {
      m_app->initBenchmark(m_argc, m_argv);
    }
    m_init = true;
  }

  int     m_argc{};
  GChar** m_argv{};

  std::unique_ptr<SolverInterface> m_app;
  GString                          m_configurationFile = "grid.json";
  GBool                            m_benchmark         = false;
  GBool                            m_init              = false;
  GInt32                           m_domainId          = -1;
  GInt32                           m_noDomains         = -1;
};
} // namespace internal_


void startupInfo(GChar** argv) {
  using namespace std;


  if(MPI::isRoot()) {
    cout << "==========================================================" << endl;
    cout << "  O)) O)   O))))))))    O))    O))       O)) O))       O))" << endl;
    cout << "O))    O)) O))       O))   O)) O) O))   O))) O) O))   O)))" << endl;
    cout << " O))       O))      O))        O)) O)) O O)) O)) O)) O O))" << endl;
    cout << "   O))     O))))))  O))        O))  O))  O)) O))  O))  O))" << endl;
    cout << "      O))  O))      O))        O))   O)  O)) O))   O)  O))" << endl;
    cout << "O))    O)) O))       O))   O)) O))       O)) O))       O))" << endl;
    cout << "  O)) O)   O))         O))))   O))       O)) O))       O))" << endl;
    cout << "                                                v" << XSTRINGIFY(PROJECT_VER) << "b" << XSTRINGIFY(BUILD_NUM) << endl;
    cout << "----------------------------------------------------------" << endl;

    cout << "Start time:            " << dateString() << "\n"
         << "Number of ranks:       " << MPI::globalNoDomains() << "\n"
#ifdef _OPENMP
         << "Number of OMP threads: " << omp_get_max_threads() << "\n"
#endif
         << "Host (of rank 0):      " << hostString() << "\n"
         << "Working directory:     " << getCWD() << "\n"
         << "Executable:            " << argv[0] << "\n"
         << endl;
  }
}

/// Default entry point of the application. Parse commandline arguments.
/// \param argc Number of received commandline arguments.
/// \param argv String of received commandline arguments.
/// \return The status of the application main run loop. (0 = ok -1 = Error...)
auto main(int argc, GChar** argv) -> int {
  std::ostringstream tmpBuffer;
#ifdef SOLVER_AVAILABLE
  tmpBuffer << "LBM Solver v" << XSTRINGIFY(PROJECT_VER);
  // todo: also give information about the version of the gridgenerator
#else
  tmpBuffer << "GridGenerator v" << XSTRINGIFY(PROJECT_VER);
#endif
  cxxopts::Options options(tmpBuffer.str(), "A highly parallel grid generator.");

  options.add_options()("d,debug", "Enable debugging with given level.", cxxopts::value<GInt>()->default_value("0"));
  options.add_options()("h,help", "Print usage");
  options.add_options()("v,version", "Get version information");
  options.add_options()("c,config", "Configuration file (default=grid.json)", cxxopts::value<std::string>()->default_value("grid.json"));
  options.add_options()("b,bench", "Run benchmark");
#ifdef SOLVER_AVAILABLE
  // todo: set solver type
  options.add_options()("s,solver", "Run solver");
#endif

  options.parse_positional({"config"});
  auto result = options.parse(argc, argv);

  if(result.count("help") > 0) {
    std::cout << options.help() << std::endl;
    exit(0);
  }
  if(result.count("version") > 0) {
    std::cout << XSTRINGIFY(PROJECT_VER) << " build " << XSTRINGIFY(BUILD_NUM) << std::endl;
    exit(0);
  }

#ifdef _OPENMP
  int provided = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
  MPI_Init(&argc, &argv);
#endif

  GInt32 domainId = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &domainId);

  GInt32 noDomains = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &noDomains);
  MPI::g_mpiInformation.init(domainId, noDomains);

  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  // Open cerr0 on MPI root
  if(MPI::isRoot()) {
    cerr0.rdbuf(std::cerr.rdbuf());
  } else {
    cerr0.rdbuf(&nullBuffer);
  }

  RESET_TIMERS();

  NEW_TIMER_GROUP_NOCREATE(TimeKeeper[Timers::AppGroup], "Application");
  NEW_TIMER_NOCREATE(TimeKeeper[Timers::timertotal], "Total", TimeKeeper[Timers::AppGroup]);
  RECORD_TIMER_START(TimeKeeper[Timers::timertotal]);

  startupInfo(argv);

  internal_::AppConfiguration gridGenRunner(argc, argv, domainId, noDomains);
#ifdef SOLVER_AVAILABLE
  internal_::AppConfiguration<SolverType::LBM> solverRunner(argc, argv, domainId, noDomains);
#endif


  GInt debug = result["debug"].as<GInt>();
  if(debug > 0) {
    if(debug > static_cast<GInt>(Debug_Level::max_debug)) {
      debug = static_cast<GInt>(Debug_Level::max_debug);
    }
    std::cout << "Activated debug level " << DEBUG_LEVEL.at(debug) << std::endl;
  }

  if(result.count("bench") > 0) {
    gridGenRunner.setBenchmark();
#ifdef SOLVER_AVAILABLE
    solverRunner.setBenchmark();
#endif
  } else {
    // first positional argument should be the configuration file
    GString config_file = result["config"].as<GString>();
    // check if the file actually exists
    if(!isFile(config_file)) {
      TERMM(-1, "Configuration file not found: " + config_file);
    }
    gridGenRunner.setConfigurationFile(config_file);
#ifdef SOLVER_AVAILABLE
    solverRunner.setConfigurationFile(config_file);
#endif
  }

  GInt ret = gridGenRunner.run(debug);
  STOP_ALL_RECORD_TIMERS();
  DISPLAY_ALL_TIMERS();
  RECORD_TIMER_START(TimeKeeper[Timers::timertotal]);


#ifdef SOLVER_AVAILABLE
  if(ret == 0 && (result.count("solver") > 0 || solverRunner.toRun(SolverType::LBM))) {
    logger.close();
    solverRunner.transferGrid(gridGenRunner.grid(), debug);
    gridGenRunner.releaseMemory();
    ret = solverRunner.run(debug);
  }
#endif

  STOP_ALL_RECORD_TIMERS();
  DISPLAY_ALL_TIMERS();

  logger.close();
  MPI_Finalize();
  return static_cast<int>(ret);
}
