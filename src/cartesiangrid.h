#ifndef GRIDGENERATOR_CARTESIANGRID_H
#define GRIDGENERATOR_CARTESIANGRID_H

#include <gcem.hpp>

#include <sfcmm_common.h>
#include "base_cartesiangrid.h"
//#include "celltree.h"
#include "common/IO.h"
#include "geometry.h"
#include "globaltimers.h"
#ifdef SOLVER_AVAILABLE
#include "gridgenerator/cartesiangrid_generation.h"
#else
#include "cartesiangrid_generation.h"
#endif
#include "interface/grid_interface.h"

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class CartesianGrid : public BaseCartesianGrid<DEBUG_LEVEL, NDIM> {
 public:
  /// Underlying enum type for property access
  using Cell = CellProperties;

  /// Underlying bitset type for property storage
  using PropertyBitsetType = grid::cell::BitsetType;

  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::checkBounds;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::property;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::size;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::empty;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::lengthOnLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::hasParent;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::parent;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::capacity;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::level;

  CartesianGrid()                     = default;
  ~CartesianGrid() override           = default;
  CartesianGrid(const CartesianGrid&) = delete;
  CartesianGrid(CartesianGrid&&)      = delete;
  auto operator=(const CartesianGrid&) -> CartesianGrid& = delete;
  auto operator=(CartesianGrid&&) -> CartesianGrid& = delete;

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
      return m_nghbrIds.at(id * cartesian::maxNoNghbrs<NDIM>() + dir);
    }
    return m_nghbrIds[id * cartesian::maxNoNghbrs<NDIM>() + dir];
  }

  [[nodiscard]] inline auto neighbor(const GInt id, const GInt dir) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
      return m_nghbrIds.at(id * cartesian::maxNoNghbrs<NDIM>() + dir);
    }
    return m_nghbrIds[id * cartesian::maxNoNghbrs<NDIM>() + dir];
  }

  [[nodiscard]] inline auto hasNeighbor(const GInt id, const GInt dir) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
      return m_nghbrIds.at(id * cartesian::maxNoNghbrs<NDIM>() + dir) > -1;
    }
    return m_nghbrIds[id * cartesian::maxNoNghbrs<NDIM>() + dir] > -1;
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

  inline auto center(const GInt id, const GInt dir) -> GDouble& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
    }
    return m_center[id][dir];
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

  /// Return number of properties defined for each node
  static constexpr auto noProperties() -> GInt { return grid::cell::p(Cell::NumProperties); }

  void setCapacity(const GInt _capacity) override {
    if(!empty()) {
      TERMM(-1, "Invalid operation tree already allocated.");
    }
    m_childIds.resize(_capacity * cartesian::maxNoChildren<NDIM>());
    m_nghbrIds.resize(_capacity * cartesian::maxNoNghbrs<NDIM>());
    m_globalIds.resize(_capacity);
    m_center.resize(_capacity);
    m_weight.resize(_capacity);
    m_noOffsprings.resize(_capacity);
    m_workload.resize(_capacity);
    BaseCartesianGrid<DEBUG_LEVEL, NDIM>::setCapacity(_capacity);
    reset();
  }

  void reset() override {
    std::fill(m_childIds.begin(), m_childIds.end(), INVALID_CELLID);
    std::fill(m_nghbrIds.begin(), m_nghbrIds.end(), INVALID_CELLID);
    std::fill(m_globalIds.begin(), m_globalIds.end(), INVALID_CELLID);
    std::fill(m_weight.begin(), m_weight.end(), NAN);
    std::fill(m_noOffsprings.begin(), m_noOffsprings.end(), INVALID_CELLID);
    std::fill(m_workload.begin(), m_workload.end(), NAN);
    for(GInt i = 0; i < capacity(); ++i) {
      property(i).reset();
      m_center[i].fill(NAN);
    }
    BaseCartesianGrid<DEBUG_LEVEL, NDIM>::reset();
  }

  void save(const GString& fileName, const json& gridOutConfig) override { TERMM(-1, "Not implemented!"); }

  /// Load the generated grid in-memory and set additional properties
  /// \param grid Generated grid.
  void loadGridInplace(const CartesianGridGen<DEBUG_LEVEL, NDIM>& grid) {
    // grid.balance(); //todo: implement
    setCapacity(grid.capacity()); // todo: change for adaptation
    m_geometry = grid.geometry();

#ifdef _OPENMP
#pragma omp parallel default(none) shared(grid)
    {
#endif
      const GInt noCells = grid.size();
#ifdef _OPENMP
#pragma omp for
#endif
      for(GInt cellId = 0; cellId < noCells; ++cellId) {
        globalId(cellId) = grid.globalId(cellId);
        parent(cellId)   = grid.parent(cellId);
        for(GInt childId = 0; childId < cartesian::maxNoChildren<NDIM>(); ++childId) {
          child(cellId, childId) = grid.child(cellId, childId);
        }
        for(GInt nghbrId = 0; nghbrId < cartesian::maxNoNghbrs<NDIM>(); ++nghbrId) {
          neighbor(cellId, nghbrId) = grid.neighbor(cellId, nghbrId);
        }
        level(cellId) = grid.level(cellId);
        for(GInt dir = 0; dir < NDIM; ++dir) {
          center(cellId, dir) = grid.center(cellId, dir);
        }
      }
#ifdef _OPENMP
    }
#endif
    setProperties();
    if(m_loadBalancing) {
      setWorkload();
      calculateOffspringsAndWeights();
    }
  }

 private:
  void setProperties() { determineBoundaryCells(); };
  void determineBoundaryCells() {
    const GInt noCells = size();
    for(GInt cellId = 0; cellId < noCells; ++cellId) {
      // is a partition cell determine for each if it can be a boundary cell (no existent parent)
      // parent has a cut with the boundary -> possible cut of child!
      if(parent(cellId) == -1 || property(parent(cellId), CellProperties::bndry)) {
        const GDouble cellLength                = lengthOnLvl(std::to_integer<GInt>(level(cellId)));
        property(cellId, CellProperties::bndry) = m_geometry->cutWithCell(center(cellId), cellLength);
      }
    }
  }
  void setWorkload() { TERMM(-1, "Not implemented!"); };
  void calculateOffspringsAndWeights() { TERMM(-1, "Not implemented!"); };

  void invalidate(const GInt begin, const GInt end) {
    std::fill(&parent(begin), &parent(end), INVALID_CELLID);
    std::fill(m_childIds.begin() + begin, m_childIds.begin() + end, INVALID_CELLID);
    std::fill(m_nghbrIds.begin() + begin, m_nghbrIds.begin() + end, INVALID_CELLID);
    std::fill(m_globalIds.begin() + begin, m_globalIds.begin() + end, INVALID_CELLID);
    std::fill(&level(0), &level(end), std::byte(-1));
    std::fill(m_center.begin() + begin, m_center.begin() + end, NAN);
    std::fill(m_weight.begin() + begin, m_weight.begin() + end, NAN);
    std::fill(m_noOffsprings.begin() + begin, m_noOffsprings.begin() + end, INVALID_CELLID);
    std::fill(m_workload.begin() + begin, m_workload.begin() + end, NAN);
    for(GInt i = 0; i < capacity(); ++i) {
      property(i).reset();
    }
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
      for(GInt j = 0; j < cartesian::maxNoNghbrs<NDIM>(); j++) {
        if(hasNeighbor(i, j)) {
          neighbor(neighbor(i, j), cartesian::oppositeDir(j)) = -1;
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
      for(GInt j = 0; j < cartesian::maxNoNghbrs<NDIM>(); ++j) {
        if(hasNeighbor(destination, j)) {
          const GInt n = neighbor(destination, j);
          if(inMovedRange(n)) {
            neighbor(destination, j) += distance;
          } else {
            neighbor(n, cartesian::oppositeDir(j)) = destination;
          }
        }
      }
    }
  }

  void checkChildPos(const GInt pos) const {
    if(pos > cartesian::maxNoChildren<NDIM>() || pos < 0) {
      TERMM(-1, "Invalid child position");
    }
  }

  void checkDir(const GInt dir) const {
    if(dir > cartesian::maxNoNghbrs<NDIM>() || dir < 0) {
      TERMM(-1, "Invalid direction");
    }
  }

  //  cartesian::Tree<DEBUG_LEVEL, NDIM> m_tree{};
  std::shared_ptr<GeometryManager<DEBUG_LEVEL, NDIM>> m_geometry;

  GBool m_loadBalancing = false;

  // Data containers
  std::vector<GInt>      m_globalIds{};
  std::vector<GInt>      m_childIds{};
  std::vector<GInt>      m_nghbrIds{};
  std::vector<GInt>      m_noOffsprings{};

  std::vector<GFloat> m_weight{};
  std::vector<GFloat> m_workload{};

  std::vector<Point<NDIM>> m_center{};
};

#endif // GRIDGENERATOR_CARTESIANGRID_H
