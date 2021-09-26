#include <cxxopts.hpp>

#include <sfcmm_common.h>
#include "app_interface.h"
#include "config.h"
#include "gridGenerator.h"
#include "lbm_solver.h"

namespace internal_ {
/// Helper object to store some general application configuration options
template <SolverType Solver = SolverType::NONE>
class AppConfiguration {
 public:
  auto run(GInt debug) -> int {
    switch(debug) {
      case static_cast<GInt>(Debug_Level::no_debug):
        return run<Debug_Level::no_debug>();
      case static_cast<GInt>(Debug_Level::min_debug):
        return run<Debug_Level::min_debug>();
      case static_cast<GInt>(Debug_Level::debug):
        return run<Debug_Level::debug>();
      case static_cast<GInt>(Debug_Level::more_debug):
        return run<Debug_Level::more_debug>();
      case static_cast<GInt>(Debug_Level::max_debug):
        [[fallthrough]];
      default:
        return run<Debug_Level::max_debug>();
    }
  }

  /// Start point of the main application.
  /// \tparam DEBUG The debug level that is activated.
  /// \return The status of the application main run loop. (0 = ok -1 = Error...)
  template <Debug_Level DEBUG>
  auto run() -> int {
    std::unique_ptr<AppInterface> App;
    switch(Solver) {
#ifdef SOLVER_AVAILABLE
      case SolverType::LBM:
        App = std::make_unique<LBMSolver<DEBUG>>();
#endif
        break;
      case SolverType::NONE:
        [[fallthrough]];
      default:
        App = std::make_unique<GridGenerator<DEBUG>>();
    }

    if(!m_benchmark) {
      App->init(m_argc, m_argv, m_configurationFile);
    } else {
      App->initBenchmark(m_argc, m_argv);
    }
    return static_cast<int>(App->run());
  }

  /// Set the commandline arguments for later processing.
  /// \param argc number of commandline arguments
  /// \param argv string of arguments
  void setCMD(int argc, GChar** argv) {
    m_argc = argc;
    m_argv = argv;
  }

  void setConfigurationFile(GString& configFile) { m_configurationFile = configFile; }
  void setBenchmark() { m_benchmark = true; }

  auto toRun(SolverType solver) const -> GBool {
    return true; // todo: implement
  }

 private:
  GChar** m_argv{};
  int     m_argc{};

  GString m_configurationFile = "grid.json";
  GBool   m_benchmark         = false;
};
} // namespace internal_


/// Default entry point of the application. Parse commandline arguments.
/// \param argc Number of received commandline arguments.
/// \param argv String of received commandline arguments.
/// \return The status of the application main run loop. (0 = ok -1 = Error...)
auto main(int argc, GChar** argv) -> int {
  std::ostringstream tmpBuffer;
#ifdef SOLVER_AVAILABLE
  tmpBuffer << "LBM Solver v" << XSTRINGIFY(PROJECT_VER);
#else
  tmpBuffer << "GridGenerator v" << XSTRINGIFY(PROJECT_VER);
#endif
  cxxopts::Options options(tmpBuffer.str(), "A highly parallel grid generator.");

  options.add_options()("d,debug", "Enable debugging with given level.", cxxopts::value<GInt>()->default_value("0"));
  options.add_options()("h,help", "Print usage");
  options.add_options()("v,version", "Get version information");
  options.add_options()("c,config", "Configuration file (default=grid.json)", cxxopts::value<std::string>()->default_value("grid.json"));
  options.add_options()("b,bench", "Run benchmark");

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

  internal_::AppConfiguration gridGenRunner{};
  gridGenRunner.setCMD(argc, argv);
#ifdef SOLVER_AVAILABLE
  internal_::AppConfiguration<SolverType::LBM> solverRunner{};
  solverRunner.setCMD(argc, argv);
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

#ifdef SOLVER_AVAILABLE
  if(ret == 0 && (result.count("solver") > 0 || solverRunner.toRun(SolverType::LBM))) {
    logger.close();
    //    solverRunner.transferGrid(gridGenRunner.grid());
    //    gridGenRunner.releaseMemory();
    ret = solverRunner.run(debug);
  }
#endif

  logger.close();
  return static_cast<int>(ret);
}
