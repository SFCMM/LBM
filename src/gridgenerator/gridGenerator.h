#ifndef GRIDGENERATOR_GRIDGENERATOR_H
#define GRIDGENERATOR_GRIDGENERATOR_H

#include <json.h>
#include <ostream>

#include <sfcmm_common.h>
#include "cartesiangrid.h"
#include "configuration.h"
#include "geometry.h"
#include "globaltimers.h"
#include "gridcell_properties.h"
#include "interface/solver_interface.h"
#include "WeightMethod.h"

using json = nlohmann::json;

template <Debug_Level DEBUG_LEVEL>
class GridGenerator : public SolverInterface, private Configuration {
 public:
  GridGenerator(GInt32 domainId, GInt32 noDomains) : m_domainId(domainId), m_noDomains(noDomains){};
  ~GridGenerator() override           = default;
  GridGenerator(const GridGenerator&) = delete;
  GridGenerator(GridGenerator&&)      = delete;
  auto operator=(const GridGenerator&) -> GridGenerator& = delete;
  auto operator=(GridGenerator&&) -> GridGenerator& = delete;


  void init(int argc, GChar** argv);
  void init(int argc, GChar** argv, GString config_file) override;
  void initBenchmark(int argc, GChar** argv) override;
  auto run() -> GInt override;
  auto grid() const -> const GridInterface& override { return *m_grid; };
  void transferGrid(const GridInterface& /*grid*/) override { TERMM(-1, "Not implemented!"); };

 private:
  int m_domainId  = -1;
  int m_noDomains = -1;

  GString m_exe;
  json                               m_geometryConfig;
  json                               m_gridOutConfig;

  void initTimers();
  void loadConfiguration();
  template <GInt nDim>
  void loadGridDefinition();
  template <GInt nDim>
  void benchmarkSetup();
  template <GInt nDim>
  void generateGrid();
  template <GInt NDIM>
  [[nodiscard]] auto inline gridGen() -> CartesianGridGen<DEBUG_LEVEL, NDIM>& {
    return *static_cast<CartesianGridGen<DEBUG_LEVEL, NDIM>*>(m_grid.get());
  }

  //  template<GInt NDIM>
  //  static std::shared_ptr<GeometryManager<DEBUG_LEVEL, NDIM>> gm(){
  //      std::shared_ptr<GeometryManager<DEBUG_LEVEL, NDIM>>(m_comm);
  //  };

  GInt m_dim        = -1;
  GInt m_maxNoCells = -1;
  // maximum size of subtree of a partitioning cell
  GInt m_maxNoOffsprings = DEFAULT_MAXNOOFFSPRINGS;
  // maximum workload of subtree of a partitioning cell
  GInt                               m_maxOffspringWorkload = DEFAULT_MAXNOOFFSPRINGS;
  GInt                               m_partitionLvl         = -1;
  GInt                               m_uniformLvl           = -1;
  GInt                               m_maxRefinementLvl     = -1;
  GBool                              m_dryRun               = false;
  GBool                              m_benchmark            = false;
  GString                            m_outputDir            = "out";
  GString                            m_outGridFilename      = "grid";
  std::unique_ptr<WeightMethod>      m_weightMethod;
  std::unique_ptr<GridInterface>     m_grid;
  std::shared_ptr<GeometryInterface> m_geometry;
};

#endif // GRIDGENERATOR_GRIDGENERATOR_H
