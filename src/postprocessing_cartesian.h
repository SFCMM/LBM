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
    if constexpr(NDIM > 1) {
      m_line = std::make_unique<Line<NDIM>>(required_config_value<NDIM>("A"), required_config_value<NDIM>("B"));

      for(GInt cellId = 0; cellId < m_cartGridData.noCells(); ++cellId) {
        if(m_cartGridData.isLeaf(cellId)) {
          const GDouble distance   = m_line->distance(m_cartGridData.center(cellId));
          const GDouble cellLength = 0.5 * m_cartGridData.cellLength(cellId);
          if(cellLength >= distance) {
            m_intersectingCells.emplace_back(cellId);
          }
        }
      }
    }

    // sort the cellIds by the coordinates
    std::sort(m_intersectingCells.begin(), m_intersectingCells.end(), [=](const GInt cellIdA, const GInt cellIdB) {
      const auto& centerA = m_cartGridData.center(cellIdA);
      const auto& centerB = m_cartGridData.center(cellIdB);
      for(GInt dir = 0; dir < NDIM; ++dir) {
        if(centerA[dir] < centerB[dir]) {
          return true;
        }
      }
      return false;
    });


    // todo: improve
    //  1. Determine if line intersects with the overall bounding box

    // 2. Find leaf cell closest to the point A

    // 3. Find all cells which intersect the line from the first intersecting leaf cell
  }

  void execute() override {
    // todo: implement
    // nothing to do since we just pass the cell list
  }

  auto output() -> std::vector<GInt>& override { return m_intersectingCells; }

 private:
  std::unique_ptr<Line<NDIM>> m_line = nullptr;

  const CartesianGridData<NDIM>& m_cartGridData{};
  std::vector<GInt>              m_intersectingCells;
};

#endif // LBM_POSTPROCESSING_CARTESIAN_H
