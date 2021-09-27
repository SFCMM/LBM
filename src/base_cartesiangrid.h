#ifndef GRIDGENERATOR_BASE_CARTESIANGRID_H
#define GRIDGENERATOR_BASE_CARTESIANGRID_H

#include <gcem.hpp>

#include <sfcmm_common.h>
#include "celltree.h"
#include "geometry.h"
#include "globaltimers.h"
#include "interface/grid_interface.h"
#include "IO.h"

struct LevelOffsetType {
 public:
  GInt begin;
  GInt end;
};

inline auto levelSize(LevelOffsetType& level) -> GInt { return level.end - level.begin; }

template <GInt NDIM>
using Point = VectorD<NDIM>;

template <GInt NDIM>
struct /*alignas(64)*/ NeighborList {
  std::array<GInt, cartesian::maxNoNghbrs<NDIM>()> n{INVALID_LIST<cartesian::maxNoNghbrs<NDIM>()>()};
};

template <GInt NDIM>
struct /*alignas(64)*/ ChildList {
  std::array<GInt, cartesian::maxNoChildren<NDIM>()> c{INVALID_LIST<cartesian::maxNoChildren<NDIM>()>()};
};

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class BaseCartesianGrid : public GridInterface {
 public:
  BaseCartesianGrid()                         = default;
  ~BaseCartesianGrid() override               = default;
  BaseCartesianGrid(const BaseCartesianGrid&) = delete;
  BaseCartesianGrid(BaseCartesianGrid&&)      = delete;
  auto operator=(const BaseCartesianGrid&) -> BaseCartesianGrid& = delete;
  auto operator=(BaseCartesianGrid&&) -> BaseCartesianGrid& = delete;

  void setBoundingBox(const BoundingBoxInterface& bbox) override {
    if(bbox.size() != NDIM) {
      TERMM(-1, "Invalid boundary box definition.");
    }

    for(GInt dir = 0; dir < NDIM; ++dir) {
      m_boundingBox.min(dir) = bbox.min(dir);
      m_boundingBox.max(dir) = bbox.max(dir);
      m_geometryExtents[dir] = gcem::abs(m_boundingBox.max(dir) - m_boundingBox.min(dir));
      // direction of largest extent will be = 0 if all extents are equal
      m_decisiveDirection    = m_geometryExtents[dir] > m_geometryExtents[m_decisiveDirection] ? dir : m_decisiveDirection;
      m_centerOfGravity[dir] = m_boundingBox.min(dir) + HALF * (m_boundingBox.max(dir) - m_boundingBox.min(dir));
    }
    m_lengthOnLevel[0] =
        (1.0 + 1.0 / gcem::pow(static_cast<GDouble>(BASE2), static_cast<GDouble>(MAX_LVL))) * m_geometryExtents[m_decisiveDirection];
    for(GInt l = 1; l < MAX_LVL; l++) {
      m_lengthOnLevel.at(l) = HALF * m_lengthOnLevel.at(l - 1);
    }
  }

  void setMaxLvl(const GInt maxLvl) override {
    logger << "set maximum grid level " << maxLvl << std::endl;
    m_maxLvl = maxLvl;
  }

  void setGeometryManager(std::shared_ptr<GeometryInterface> geom) override {
    // we cast this here since we don't want to cast all the time to get rid of NDIM...
    m_geometry = std::static_pointer_cast<GeometryManager<DEBUG_LEVEL, NDIM>>(geom);
  }

  [[nodiscard]] inline auto cog() const -> std::vector<GDouble> override {
    return std::vector<GDouble>(m_centerOfGravity.begin(), m_centerOfGravity.end());
  };
  [[nodiscard]] inline auto lengthOfBoundingBox() const -> std::vector<GDouble> override {
    return std::vector<GDouble>(m_geometryExtents.begin(), m_geometryExtents.end());
  };
  [[nodiscard]] inline auto boundingBox() const -> const BoundingBoxInterface& override { return m_boundingBox; };
  [[nodiscard]] inline auto largestDir() const -> GInt override { return m_decisiveDirection; };
  [[nodiscard]] inline auto partitionLvl() const -> GInt override { return m_partitioningLvl; };
  inline auto               partitionLvl() -> GInt& override { return m_partitioningLvl; }
  [[nodiscard]] inline auto maxLvl() const -> GInt override { return m_maxLvl; };
  [[nodiscard]] inline auto lengthOnLvl(const GInt lvl) const -> GDouble override {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      return m_lengthOnLevel.at(lvl);
    }
    return m_lengthOnLevel[lvl]; // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
  };
  [[nodiscard]] inline auto currentHighestLvl() const -> GInt override { return m_currentHighestLvl; }

 protected:
  /// Increase the current highest level by 1
  inline void increaseCurrentHighestLvl() {
    ASSERT(m_currentHighestLvl <= m_maxLvl, "Level increased over maximum level!");
    ++m_currentHighestLvl;
  }

  /// Give write access to geometry
  inline auto geometry() { return m_geometry; }

  /// Get access to geometry
  inline auto geometry() const { return m_geometry; }

 private:
  std::shared_ptr<GeometryManager<DEBUG_LEVEL, NDIM>> m_geometry;

  GInt m_currentHighestLvl = 0;
  GInt m_partitioningLvl   = 0;
  GInt m_maxLvl            = 0;

  // box containing the whole geometry
  BoundingBoxCT<NDIM> m_boundingBox;
  //  std::array<GDouble, 2 * NDIM> m_boundingBox{NAN_LIST<2 * NDIM>()};
  // extent of the geometry
  std::array<GDouble, NDIM> m_geometryExtents{NAN_LIST<NDIM>()};
  // m_center of gravity of the geometry
  std::array<GDouble, NDIM> m_centerOfGravity{NAN_LIST<NDIM>()};
  // direction of largest extent
  GInt m_decisiveDirection{};
  // length of the cells on each level basest on the largest extent
  std::array<GDouble, MAX_LVL> m_lengthOnLevel{NAN_LIST<MAX_LVL>()};
};

#endif // GRIDGENERATOR_CARTESIANGRID_H
