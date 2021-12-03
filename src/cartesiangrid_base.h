#ifndef GRIDGENERATOR_BASE_CARTESIANGRID_H
#define GRIDGENERATOR_BASE_CARTESIANGRID_H

#include <gcem.hpp>

#include <sfcmm_common.h>
//#include "celltree.h"
#include "common/IO.h"
#include "geometry.h"
#include "globaltimers.h"
#include "gridcell_properties.h"
#include "interface/grid_interface.h"

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
class BaseCartesianGrid;

/// Data access object to give const-level access to an existing cartesian grid.
/// \tparam NDIM
template <GInt NDIM>
class CartesianGridData {
 private:
  using PropertyBitsetType = grid::cell::BitsetType;

 public:
  template <typename T>
  CartesianGridData(const T& initGrid)
    : m_noCells(initGrid.noCells()),
      m_boundingBox(initGrid.boundingBox()),
      m_center(initGrid.center()),
      m_properties(initGrid.props()),
      m_level(initGrid.level()),
      m_lengthOnLevel(initGrid.lengthOnLvl()) {}

  // todo: add asserts
  [[nodiscard]] inline auto noCells() const -> GInt { return m_noCells; }

  // todo: add asserts
  inline auto center(const GInt cellId) const -> const Point<NDIM>& { return m_center[cellId]; }

  // todo: add asserts
  [[nodiscard]] inline auto isLeaf(const GInt cellId) const -> GBool {
    return m_properties[cellId][static_cast<GInt>(CellProperties::leaf)];
  }

  // todo: add asserts
  [[nodiscard]] inline auto level(const GInt cellId) const -> std::byte { return m_level[cellId]; }

  // todo: add asserts
  [[nodiscard]] inline auto cellLength(const GInt cellId) const -> GDouble { return m_lengthOnLevel[static_cast<GInt>(level(cellId))]; }

  [[nodiscard]] inline auto boundingBox() const -> const BoundingBoxInterface& { return m_boundingBox; }


 private:
  const GInt m_noCells;

  const BoundingBoxInterface&            m_boundingBox;
  const std::vector<Point<NDIM>>&        m_center;
  const std::vector<PropertyBitsetType>& m_properties;
  const std::vector<std::byte>&          m_level;
  const std::array<GDouble, MAX_LVL>     m_lengthOnLevel{NAN_LIST<MAX_LVL>()};
};

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class BaseCartesianGrid : public GridInterface {
 public:
  using PropertyBitsetType = grid::cell::BitsetType;

  /// Underlying enum type for property access
  using Cell = CellProperties;

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
  [[nodiscard]] inline auto lengthOnLvl(const std::byte lvl) const -> GDouble override { return lengthOnLvl(static_cast<GInt>(lvl)); };
  [[nodiscard]] inline auto lengthOnLvl() const -> const auto& { return m_lengthOnLevel; };

  [[nodiscard]] inline auto cellLength(const GInt cellId) const -> GDouble override {
    const GInt lvl = static_cast<GInt>(level(cellId));
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      return m_lengthOnLevel.at(lvl);
    }
    return m_lengthOnLevel[lvl]; // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
  };
  [[nodiscard]] inline auto currentHighestLvl() const -> GInt override { return m_currentHighestLvl; }
  [[nodiscard]] inline auto dim() const -> GInt override { return NDIM; }

  inline auto property(const GInt id, CellProperties p) -> auto {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      //      checkBounds(id);
      checkProperty(p);
      return m_properties.at(id)[grid::cell::p(p)];
    }
    return m_properties[id][static_cast<GInt>(p)];
  }
  [[nodiscard]] inline auto property(const GInt id, CellProperties p) const -> GBool override {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkProperty(p);
      return m_properties.at(id)[grid::cell::p(p)];
    }
    return m_properties[id][static_cast<GInt>(p)];
  }

  [[nodiscard]] inline auto center(const GInt id, const GInt dir) const -> GDouble {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
    }
    return m_center[id][dir];
  }

  [[nodiscard]] inline auto center(const GInt id) const -> const Point<NDIM>& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
    }
    return m_center[id];
  }

  inline auto center() const -> const std::vector<Point<NDIM>>& { return m_center; }


  [[nodiscard]] inline auto level(const GInt id) const -> std::byte override {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_level.at(id);
    }
    return m_level[id];
  }

  [[nodiscard]] inline auto level() -> std::vector<std::byte>& { return m_level; }
  [[nodiscard]] inline auto level() const -> const std::vector<std::byte>& { return m_level; }

  [[nodiscard]] inline auto globalId(const GInt id) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_globalId.at(id);
    }
    return m_globalId[id];
  }

  [[nodiscard]] inline auto propertiesToBits(const GInt id) const -> GUint {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_properties.at(id).to_ullong();
    }
    return m_properties[id].to_ullong();
  }

  void propertiesFromBits(const GInt id, const GUint bits) {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
    }
    m_properties[id] = PropertyBitsetType(bits);
  }

  [[nodiscard]] inline auto propertiesToString(const GInt id) const -> GString {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_properties.at(id).to_string();
    }
    return m_properties[id].to_string();
  }

  inline auto property(const GInt id) -> PropertyBitsetType& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      //      checkBounds(id);
      return m_properties.at(id);
    }
    return m_properties[id];
  }

  void resetProperties(const GInt id) {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
    }
    m_properties[id].reset();
  }

  [[nodiscard]] inline auto isLeafCell(const GInt id) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_properties.at(id)[static_cast<GInt>(CellProperties::leaf)];
    }
    return m_properties[id][static_cast<GInt>(CellProperties::leaf)];
  }

  void checkBounds(const GInt id) const {
    if(id > this->capacity()) {
      TERMM(-1, "Out of bounds id: " + std::to_string(id) + "/" + std::to_string(this->capacity()));
    }
  }

  [[nodiscard]] inline auto capacity() const -> GInt { return m_capacity; }
  [[nodiscard]] inline auto size() const -> GInt override { return m_size; }
  [[nodiscard]] inline auto noCells() const -> GInt override { return m_size; }
  [[nodiscard]] inline auto empty() const -> GBool { return m_size == 0; }

  [[nodiscard]] inline auto parent(const GInt id) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_parentId.at(id);
    }
    // no bound checking
    return m_parentId[id];
  }

  [[nodiscard]] inline auto hasParent(const GInt id) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_parentId.at(id) > -1;
    }
    // no bound checking
    return m_parentId[id] > -1;
  }

  [[nodiscard]] inline auto props() const -> const std::vector<PropertyBitsetType>& { return m_properties; }


  auto getCartesianGridData() const -> CartesianGridData<NDIM> { return CartesianGridData<NDIM>(*this); }

 protected:
  inline auto currentHighestLvl() -> GInt& { return m_currentHighestLvl; }

  inline auto parent(const GInt id) -> GInt& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      //      checkBounds(id);
      return m_parentId.at(id);
    }
    // no bound checking
    return m_parentId[id];
  }

  inline auto center() -> std::vector<Point<NDIM>>& { return m_center; }

  [[nodiscard]] inline auto center(const GInt id) -> Point<NDIM>& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      m_center.at(id);
    }
    return m_center[id];
  }

  inline auto center(const GInt id, const GInt dir) -> GDouble& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      //      checkBounds(id);
      checkDir(dir);
    }
    return m_center[id][dir];
  }

  inline auto globalId(const GInt id) -> GInt& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      //      checkBounds(id);
      return m_globalId.at(id);
    }
    return m_globalId[id];
  }

  inline auto level(const GInt id) -> std::byte& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      //      checkBounds(id);
      return m_level.at(id);
    }
    return m_level[id];
  }

  inline auto size() -> GInt& { return m_size; }

  /// Increase the current highest level by 1
  inline void increaseCurrentHighestLvl() {
    ASSERT(m_currentHighestLvl <= m_maxLvl, "Level increased over maximum level!");
    ++m_currentHighestLvl;
  }

  /// Give write access to geometry
  inline auto geometry() { return m_geometry; }

  /// Get access to geometry
  inline auto geometry() const { return m_geometry; }

  void setCapacity(const GInt capacity) override {
    m_properties.resize(capacity);
    m_center.resize(capacity);
    m_parentId.resize(capacity);
    m_level.resize(capacity);
    m_globalId.resize(capacity);
    m_capacity = capacity;
  }

  void reset() override {
    std::for_each(m_properties.begin(), m_properties.end(), [](auto& prop) { prop.reset(); });
    std::fill(m_parentId.begin(), m_parentId.end(), INVALID_CELLID);
    std::fill(m_level.begin(), m_level.end(), std::byte(-1));
    std::fill(m_globalId.begin(), m_globalId.end(), INVALID_CELLID);
    m_size = 0;
  }

  void clear() {
    m_properties.clear();
    m_parentId.clear();
    m_level.clear();
    m_center.clear();
    m_globalId.clear();
    m_size     = 0;
    m_capacity = 0;
    m_size     = 0;
  }

  void checkDir(const GInt dir) const {
    if(dir > cartesian::maxNoNghbrs<NDIM>() || dir < 0) {
      TERMM(-1, "Invalid direction " + std::to_string(dir));
    }
  }


 private:
  void checkProperty(const Cell p) const {
    if(p >= Cell::NumProperties) {
      TERMM(-1, "Invalid property!");
    }
  }

  std::shared_ptr<GeometryManager<DEBUG_LEVEL, NDIM>> m_geometry;

  GInt m_currentHighestLvl = 0;
  GInt m_partitioningLvl   = 0;
  GInt m_maxLvl            = 0;
  GInt m_size              = 0;
  GInt m_capacity          = 0;


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

  std::vector<PropertyBitsetType> m_properties{};
  std::vector<GInt>               m_parentId{};
  std::vector<GInt>               m_globalId{};
  std::vector<Point<NDIM>>        m_center{};
  std::vector<std::byte>          m_level{};
};

#endif // GRIDGENERATOR_BASE_CARTESIANGRID_H
