#ifndef LBM_POSTPROCESSING_CARTESIAN_H
#define LBM_POSTPROCESSING_CARTESIAN_H
#include <common/mesh.h>
#include <utility>
#include "common/line.h"
#include "postprocessing_func.h"

template <GInt NDIM>
class PostprocessCartesianFunctionLine : public PostprocessFunctionInterface<NDIM>, Configuration {
 public:
  PostprocessCartesianFunctionLine(json conf, const CartesianGridData<NDIM>& data) : m_cartGridData(data) {
    Configuration::setConfiguration(std::move(conf));

    init();
  };

  void init() override {
    m_line = std::make_unique<Line<NDIM>>(required_config_value<NDIM>("A"), required_config_value<NDIM>("B"));

    for(GInt cellId = 0; cellId < m_cartGridData.noCells(); ++cellId) {
      cerr0 << m_cartGridData.center(cellId) << std::endl;
      cerr0 << m_line->distance(m_cartGridData.center(cellId)) << std::endl;
    }


    // todo: improve
    //  1. Determine if line intersects with the overall bounding box

    // 2. Find leaf cell closest to the point A

    // 3. Find all cells which intersect the line from the first intersecting leaf cell
  }

  void execute() override {}

 private:
  std::unique_ptr<Line<NDIM>> m_line = nullptr;

  const CartesianGridData<NDIM>&   m_cartGridData{};
  std::vector<GInt> m_intersectingCells;
};

#endif // LBM_POSTPROCESSING_CARTESIAN_H
