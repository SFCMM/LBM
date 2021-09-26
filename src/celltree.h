#ifndef GRIDGENERATOR_CELLTREE_H
#define GRIDGENERATOR_CELLTREE_H

#include <bitset>
#include <sfcmm_common.h>
#include "functions.h"
#include "gridcell_properties.h"


namespace cartesian {
/// Underlying enum type for property access
using Cell = CellProperties;

/// Underlying bitset type for property storage
using PropertyBitsetType = grid::cell::BitsetType;

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class Tree {
 public:
  constexpr Tree() = default;

  // Parent-child relationship
  inline auto parent(const GInt id) -> GInt& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_parentIds.at(id);
    }
    // no bound checking
    return m_parentIds[id];
  }

  [[nodiscard]] inline auto parent(const GInt id) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_parentIds.at(id);
    }
    // no bound checking
    return m_parentIds[id];
  }

  [[nodiscard]] inline auto hasParent(const GInt id) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_parentIds.at(id) > -1;
    }
    // no bound checking
    return m_parentIds[id] > -1;
  }

  inline auto child(const GInt id, const GInt pos) -> GInt& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkChildPos(pos);
      checkBounds(id);
      return m_childIds.at(id * cartesian::maxNoChildren<NDIM>() + pos);
    }
    // no bound checking
    return m_childIds[id * cartesian::maxNoChildren<NDIM>() + pos];
  }

  [[nodiscard]] inline auto child(const GInt id, const GInt pos) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkChildPos(pos);
      checkBounds(id);
      return m_childIds.at(id * cartesian::maxNoChildren<NDIM>() + pos);
    }
    // no bound checking
    return m_childIds[id * cartesian::maxNoChildren<NDIM>() + pos];
  }

  [[nodiscard]] inline auto hasChild(const GInt id, const GInt pos) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkChildPos(pos);
      checkBounds(id);
      return m_childIds.at(id * cartesian::maxNoChildren<NDIM>() + pos) > -1;
    }
    // no bound checking
    return m_childIds[id * cartesian::maxNoChildren<NDIM>() + pos] > -1;
  }

  [[nodiscard]] inline auto hasChildren(const GInt id) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
    }
    return noChildren(id) > 0;
  }

  [[nodiscard]] inline auto noChildren(const GInt id) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
    }
    return std::count_if(&m_childIds[id * cartesian::maxNoChildren<NDIM>() + 0],
                         &m_childIds[id * cartesian::maxNoChildren<NDIM>() + cartesian::maxNoChildren<NDIM>()],
                         [](const GInt childId) { return childId > -1; });
  }

  // Neighbors
  inline auto neighbor(const GInt id, const GInt dir) -> GInt& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
      return m_neighborIds.at(id * noNeighborsPerNode() + dir);
    }
    return m_neighborIds[id * noNeighborsPerNode() + dir];
  }

  [[nodiscard]] inline auto neighbor(const GInt id, const GInt dir) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
      return m_neighborIds.at(id * noNeighborsPerNode() + dir);
    }
    return m_neighborIds[id * noNeighborsPerNode() + dir];
  }

  [[nodiscard]] inline auto hasNeighbor(const GInt id, const GInt dir) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
      return m_neighborIds.at(id * noNeighborsPerNode() + dir) > -1;
    }
    return m_neighborIds[id * noNeighborsPerNode() + dir] > -1;
  }

  [[nodiscard]] inline auto hasAnyNeighbor(const GInt id, const GInt dir) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
    }
    return hasNeighbor(id, dir) || (hasParent(id) && hasNeighbor(parent(id), dir));
  }

  // Other data fields
  inline auto globalId(const GInt id) -> GInt& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_globalIds.at(id);
    }
    return m_globalIds[id];
  }

  [[nodiscard]] inline auto globalId(const GInt id) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_globalIds.at(id);
    }
    return m_globalIds[id];
  }

  inline auto level(const GInt id) -> GInt& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_levels.at(id);
    }
    return m_levels[id];
  }

  [[nodiscard]] inline auto level(const GInt id) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_levels.at(id);
    }
    return m_levels[id];
  }

  inline auto coordinate(const GInt id, const GInt dir) -> GDouble& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
      return m_coordinates.at(id * NDIM + dir);
    }
    return m_coordinates[id * NDIM + dir];
  }

  [[nodiscard]] inline auto coordinate(const GInt id, const GInt dir) const -> GDouble {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
      return m_coordinates.at(id * NDIM + dir);
    }
    return m_coordinates[id * NDIM + dir];
  }

  [[nodiscard]] inline auto coordinate(const GInt id) const -> const GDouble* {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return &m_coordinates.at(id * NDIM);
    }
    return &m_coordinates[id * NDIM];
  }

  inline auto weight(const GInt id) -> GFloat& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_weight.at(id);
    }
    return m_weight[id];
  }

  [[nodiscard]] inline auto weight(const GInt id) const -> GFloat {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_weight.at(id);
    }
    return m_weight[id];
  }


  [[nodiscard]] inline auto isLeafCell(const GInt id) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_properties.at(id)[static_cast<GInt>(CellProperties::leaf)];
    }
    return m_properties[id][static_cast<GInt>(CellProperties::leaf)];
  }

  inline auto hasProperty(const GInt id, const Cell p) -> PropertyBitsetType::reference {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkProperty(p);
      return m_properties.at(id)[grid::cell::p(p)];
    }
    return m_properties[id][grid::cell::p(p)];
  }

  [[nodiscard]] inline auto hasProperty(const GInt id, const Cell p) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkProperty(p);
      return m_properties.at(id)[grid::cell::p(p)];
    }
    return m_properties[id][grid::cell::p(p)];
  }

  void resetProperties(const GInt id) {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
    }
    m_properties[id].reset();
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
  inline auto properties(const GInt id) -> PropertyBitsetType& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_properties.at(id);
    }
    return m_properties[id];
  }

  // Other data fields (subject to change)
  inline auto noOffsprings(const GInt id) -> GInt& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_noOffsprings.at(id);
    }
    return m_noOffsprings[id];
  }

  [[nodiscard]] inline auto noOffsprings(const GInt id) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_noOffsprings.at(id);
    }
    return m_noOffsprings[id];
  }

  inline auto workload(const GInt id) -> GFloat& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_workload.at(id);
    }
    return m_workload[id];
  }

  [[nodiscard]] inline auto workload(const GInt id) const -> GFloat {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      return m_workload.at(id);
    }
    return m_workload[id];
  }

  // Entries per tree node
  /// Return maximum number of same-level neighbors per node
  static constexpr auto noNeighborsPerNode() -> GInt { return 2 * NDIM; }

  /// Return opposite direction for given neighbor direction
  static constexpr auto oppositeNeighborDir(const GInt dir) -> GInt { return dir + 1 - 2 * (dir % 2); }

  /// Return number of properties defined for each node
  static constexpr auto noProperties() -> GInt { return grid::cell::p(Cell::NumProperties); }

  [[nodiscard]] inline auto capacity() const -> GInt { return m_parentIds.capacity(); }
  [[nodiscard]] inline auto size() const -> GInt { return m_size; }
  [[nodiscard]] inline auto empty() const -> GBool { return m_size == 0; }

  void reset(const GInt capacity) {
    m_parentIds.resize(capacity);
    m_childIds.resize(capacity);
    m_neighborIds.resize(capacity);
    m_globalIds.resize(capacity);
    m_levels.resize(capacity);
    m_coordinates.resize(capacity);
    m_weight.resize(capacity);
    m_properties.resize(capacity);
    m_noOffsprings.resize(capacity);
    m_workload.resize(capacity);
    reset();
  }

 private:
  void reset() {
    std::fill(m_parentIds.begin(), m_parentIds.end(), -1);
    std::fill(m_childIds.begin(), m_childIds.end(), -1);
    std::fill(m_neighborIds.begin(), m_neighborIds.end(), -1);
    std::fill(m_globalIds.begin(), m_globalIds.end(), -1);
    std::fill(m_levels.begin(), m_levels.end(), -1);
    std::fill(m_coordinates.begin(), m_coordinates.end(), NAN);
    std::fill(m_weight.begin(), m_weight.end(), NAN);
    std::fill(m_properties.begin(), m_properties.end(), 0);
    std::fill(m_noOffsprings.begin(), m_noOffsprings.end(), -1);
    std::fill(m_workload.begin(), m_workload.end(), NAN);
    m_size = 0;
  }

  void invalidate(const GInt begin, const GInt end) {
    std::fill(m_parentIds.begin() + begin, m_parentIds.begin() + end, -1);
    std::fill(m_childIds.begin() + begin, m_childIds.begin() + end, -1);
    std::fill(m_neighborIds.begin() + begin, m_neighborIds.begin() + end, -1);
    std::fill(m_globalIds.begin() + begin, m_globalIds.begin() + end, -1);
    std::fill(m_levels.begin() + begin, m_levels.begin() + end, -1);
    std::fill(m_coordinates.begin() + begin, m_coordinates.begin() + end, NAN);
    std::fill(m_weight.begin() + begin, m_weight.begin() + end, NAN);
    std::fill(m_properties.begin() + begin, m_properties.begin() + end, 0);
    std::fill(m_noOffsprings.begin() + begin, m_noOffsprings.begin() + end, -1);
    std::fill(m_workload.begin() + begin, m_workload.begin() + end, NAN);
  }
  //  template <class Functor, class T>
  //  void rawCopyGeneric(Functor&& c, const T& source, const GInt begin, const GInt end, const GInt destination);
  void deleteConnectivity(const GInt begin, const GInt end) {
    for(GInt i = begin; i < end; i++) {
      // Parent
      if(hasParent(i)) {
        const GInt p = parent(i);
        for(GInt j = 0; j < cartesian::maxNoChildren<NDIM>(); j++) {
          if(child(p, j) == i) {
            child(p, j) = -1;
          }
        }
      }

      // Children
      for(GInt j = 0; j < cartesian::maxNoChildren<NDIM>(); j++) {
        if(hasChild(i, j)) {
          parent(child(i, j)) = -1;
        }
      }

      // Neighbors
      for(GInt j = 0; j < noNeighborsPerNode(); j++) {
        if(hasNeighbor(i, j)) {
          neighbor(neighbor(i, j), oppositeNeighborDir(j)) = -1;
        }
      }
    }
  }
  void moveConnectivity(const GInt begin, const GInt end, const GInt to) {
    // Auxiliary method for checking if a given id is within the original range that was moved
    auto inMovedRange = [begin, end](const GInt id) { return (id >= begin && id < end); };

    // General strategy:
    // 1) Loop over moved nodes and check all tree connections (parents/children/neighbors)
    // 2) If a given connection is to a node that was moved: apply offset to current node
    // 3) If a given connection is to a node that was not moved: change connectivity in other node
    for(GInt from = begin; from < end; ++from) {
      const GInt distance    = to - begin;
      const GInt destination = from + distance;

      // Parent
      if(hasParent(destination)) {
        const GInt p = parent(destination);
        if(inMovedRange(p)) {
          parent(destination) += distance;
        } else {
          for(GInt j = 0; j < cartesian::maxNoChildren<NDIM>(); ++j) {
            if(child(p, j) == from) {
              child(p, j) = destination;
            }
          }
        }
      }

      // Children
      for(GInt j = 0; j < cartesian::maxNoChildren<NDIM>(); ++j) {
        if(hasChild(destination, j)) {
          const GInt c = child(destination, j);
          if(inMovedRange(c)) {
            child(destination, j) += distance;
          } else {
            parent(c) = destination;
          }
        }
      }

      // Neighbors
      for(GInt j = 0; j < noNeighborsPerNode(); ++j) {
        if(hasNeighbor(destination, j)) {
          const GInt n = neighbor(destination, j);
          if(inMovedRange(n)) {
            neighbor(destination, j) += distance;
          } else {
            neighbor(n, oppositeNeighborDir(j)) = destination;
          }
        }
      }
    }
  }

  void checkBounds(const GInt id) const {
    if(id > size()) {
      TERMM(-1, "Out of bounds.");
    }
  }

  void checkChildPos(const GInt pos) const {
    if(pos > cartesian::maxNoChildren<NDIM>() || pos < 0) {
      TERMM(-1, "Invalid child position");
    }
  }

  void checkDir(const GInt dir) const {
    if(dir > noNeighborsPerNode() || dir < 0) {
      TERMM(-1, "Invalid direction");
    }
  }

  void checkProperty(const Cell p) const {
    if(p != Cell::NumProperties) {
      TERMM(-1, "Invalid property!");
    }
  }

  // Data containers
  std::vector<GInt> m_globalIds{};
  std::vector<GInt> m_parentIds{};
  std::vector<GInt> m_childIds{};
  std::vector<GInt> m_neighborIds{};
  std::vector<GInt> m_levels{};
  std::vector<GInt> m_noOffsprings{};

  std::vector<PropertyBitsetType> m_properties{};

  std::vector<GFloat> m_weight{};
  std::vector<GFloat> m_workload{};

  std::vector<GDouble> m_coordinates{};

  GInt m_size = 0;
};
} // namespace cartesian

#endif // GRIDGENERATOR_CELLTREE_H
