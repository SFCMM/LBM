#ifndef GRIDGENERATOR_CARTESIANGRID_H
#define GRIDGENERATOR_CARTESIANGRID_H

#include <gcem.hpp>

#include <common/surface.h>
#include <set>
#include <sfcmm_common.h>
#include "cartesiangrid_base.h"
#include "common/configuration.h"
#include "common/IO.h"
#include "geometry.h"
#include "globaltimers.h"
#ifdef SOLVER_AVAILABLE
#include "gridgenerator/cartesiangrid_generation.h"
#include "lbm/lbm_constants.h"
#else
#include "cartesiangrid_generation.h"
#endif
#include "interface/grid_interface.h"

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class CartesianGrid : public BaseCartesianGrid<DEBUG_LEVEL, NDIM>, private Configuration {
 public:
  /// Underlying enum type for property access
  using Cell = CellProperties;

  /// Underlying bitset type for property storage
  using PropertyBitsetType = grid::cell::BitsetType;

  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::checkBounds;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::property;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::size;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::noCells;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::empty;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::lengthOnLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::hasParent;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::parent;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::capacity;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::level;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::globalId;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::center;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::checkDir;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::currentHighestLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::partitionLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::setMaxLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::boundingBox;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::setBoundingBox;

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
  inline auto neighbor() const -> const auto& { return m_nghbrIds; }

  // todo: make diagonal neighbors a template parameter
  [[nodiscard]] inline auto neighbor(const GInt id, const GInt dir) const -> GInt override {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      //      checkDir(dir);
      if(m_diagonalNghbrs) {
        return m_nghbrIds.at(id * cartesian::maxNoNghbrsDiag<NDIM>() + dir);
      }
      return m_nghbrIds.at(id * cartesian::maxNoNghbrs<NDIM>() + dir);
    }
    if(m_diagonalNghbrs) {
      return m_nghbrIds[id * cartesian::maxNoNghbrsDiag<NDIM>() + dir];
    }
    return m_nghbrIds[id * cartesian::maxNoNghbrs<NDIM>() + dir];
  }

  inline auto neighbor(const GInt id, const GInt dir) -> GInt& {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      //      checkDir(dir);
      if(m_diagonalNghbrs) {
        return m_nghbrIds.at(id * cartesian::maxNoNghbrsDiag<NDIM>() + dir);
      }
      return m_nghbrIds.at(id * cartesian::maxNoNghbrs<NDIM>() + dir);
    }
    if(m_diagonalNghbrs) {
      return m_nghbrIds[id * cartesian::maxNoNghbrsDiag<NDIM>() + dir];
    }
    return m_nghbrIds[id * cartesian::maxNoNghbrs<NDIM>() + dir];
  }


  [[nodiscard]] inline auto hasNeighbor(const GInt id, const GInt dir) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
      if(m_diagonalNghbrs) {
        return m_nghbrIds.at(id * cartesian::maxNoNghbrsDiag<NDIM>() + dir) != INVALID_CELLID;
      }
      return m_nghbrIds.at(id * cartesian::maxNoNghbrs<NDIM>() + dir) != INVALID_CELLID;
    }
    if(m_diagonalNghbrs) {
      return m_nghbrIds[id * cartesian::maxNoNghbrsDiag<NDIM>() + dir] != INVALID_CELLID;
    }
    return m_nghbrIds[id * cartesian::maxNoNghbrs<NDIM>() + dir] != INVALID_CELLID;
  }

  [[nodiscard]] inline auto hasAnyNeighbor(const GInt id, const GInt dir) const -> GBool {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      checkBounds(id);
      checkDir(dir);
    }
    return hasNeighbor(id, dir) || (hasParent(id) && hasNeighbor(parent(id), dir));
  }

  // Other data fields
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

  [[nodiscard]] inline auto noLeafCells() const -> GInt { return m_noLeafCells; }

  [[nodiscard]] inline auto noBndCells() const -> GInt { return m_noBndCells; }

  void setCapacity(const GInt _capacity) override {
    if(!empty()) {
      TERMM(-1, "Invalid operation tree already allocated.");
    }
    m_childIds.resize(_capacity * cartesian::maxNoChildren<NDIM>());
    m_nghbrIds.resize(_capacity * cartesian::maxNoNghbrsDiag<NDIM>());
    m_weight.resize(_capacity);
    m_noOffsprings.resize(_capacity);
    m_workload.resize(_capacity);
    BaseCartesianGrid<DEBUG_LEVEL, NDIM>::setCapacity(_capacity);
    reset();
  }

  void reset() override {
    std::fill(m_childIds.begin(), m_childIds.end(), INVALID_CELLID);
    std::fill(m_nghbrIds.begin(), m_nghbrIds.end(), INVALID_CELLID);
    std::fill(m_weight.begin(), m_weight.end(), NAN);
    std::fill(m_noOffsprings.begin(), m_noOffsprings.end(), INVALID_CELLID);
    std::fill(m_workload.begin(), m_workload.end(), NAN);
    for(GInt i = 0; i < capacity(); ++i) {
      property(i).reset();
      center(i).fill(NAN);
    }
    BaseCartesianGrid<DEBUG_LEVEL, NDIM>::reset();
  }

  void save(const GString& /*fileName*/, const json& /*gridOutConfig*/) const override { TERMM(-1, "Not implemented!"); }

  auto bndrySurface(const GString& id) -> Surface<DEBUG_LEVEL, NDIM>& {
    if(DEBUG_LEVEL > Debug_Level::min_debug) {
      if(m_bndrySurfaces.count(id) == 0) {
        TERMM(-1, "Invalid bndryId \"" + id + "\"");
      }
    }
    return m_bndrySurfaces.at(id);
  }


  /// Load the generated grid in-memory and set additional properties
  /// \param grid Generated grid.
  void loadGridInplace(const CartesianGridGen<DEBUG_LEVEL, NDIM>& grid, const json& properties) {
    setConfiguration(properties);
    // grid.balance(); //todo: implement
    setCapacity(grid.capacity()); // todo: change for adaptation
    m_geometry          = grid.geometry();
    size()              = grid.size();
    currentHighestLvl() = grid.currentHighestLvl();
    partitionLvl()      = grid.partitionLvl();
    setMaxLvl(grid.maxLvl());
    setBoundingBox(grid.boundingBox());

#ifdef _OPENMP
#pragma omp parallel default(none) shared(grid)
    {
#endif
#ifdef _OPENMP
#pragma omp for
#endif
      for(GInt cellId = 0; cellId < noCells(); ++cellId) {
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

    m_axisAlignedBnd = opt_config_value<GBool>("assumeAxisAligned", m_axisAlignedBnd);
    if(!m_axisAlignedBnd) {
      TERMM(-1, "Not implemented!");
    }
    m_periodic = has_any_key_value("type", "periodic");


    setProperties();

    determineBoundaryCells();
    identifyBndrySurfaces();
    setupPeriodicConnections();
    addGhostCells();
    addDiagonalNghbrs();
    //    if(m_loadBalancing) {
    //      setWorkload();
    //      calculateOffspringsAndWeights();
    //    }
  }

  /// Add ghost cells
  void addGhostCells() {
    // todo: make settable
    const GBool addGhostLayers = false;
    if(addGhostLayers) {
      // check all surfaces and add ghostcells in all missing dist directions
      for(const auto& [srfName, srf] : m_bndrySurfaces) {
        for(GInt cellId : srf.getCellList()) {
          for(GInt nghbrDir = 0; nghbrDir < cartesian::maxNoNghbrs<NDIM>(); ++nghbrDir) {
            if(neighbor(cellId, nghbrDir) == INVALID_CELLID) {
              const GInt ghostCellId = size() + m_noGhostsCells;
              cerr0 << "adding cell " << ghostCellId << " as neighbor to " << cellId << std::endl; // todo: remove
              neighbor(cellId, nghbrDir)                              = ghostCellId;
              neighbor(ghostCellId, cartesian::oppositeDir(nghbrDir)) = cellId;
              property(ghostCellId, CellProperties::ghost);
              ++m_noGhostsCells;
            }
          }
        }
      }
    }
  }

  /// Add the diagonal(2D/3D) and/or tridiagonal (3D) to the neighbor connections of each cell.
  void addDiagonalNghbrs() {
    m_diagonalNghbrs = true;
    auto tmpNghbr    = m_nghbrIds;

    auto tmpN = [&](const GInt cellId, const GInt dir) { return tmpNghbr[cellId * cartesian::maxNoNghbrs<NDIM>() + dir]; };


    for(GInt cellId = 0; cellId < size(); ++cellId) {
      // dirs 0=-x 1=+x 2=-y 3=+y

      // copy existing neighbor connections
      for(GInt dir = 0; dir < cartesian::maxNoNghbrs<NDIM>(); ++dir) {
        neighbor(cellId, dir) = tmpN(cellId, dir);
      }

      //  add diagonal nghbrs
      if(NDIM > 1) {
        const GInt nghbrmX = tmpN(cellId, 0);
        const GInt nghbrpX = tmpN(cellId, 1);

        // +x+y
        const GInt nghbrpXpY                             = nghbrpX != INVALID_CELLID ? tmpN(nghbrpX, 3) : -1;
        neighbor(cellId, cartesian::maxNoNghbrs<NDIM>()) = nghbrpXpY;

        // +x-y
        const GInt nghbrpXmY                                 = nghbrpX != INVALID_CELLID ? tmpN(nghbrpX, 2) : -1;
        neighbor(cellId, cartesian::maxNoNghbrs<NDIM>() + 1) = nghbrpXmY;

        // -x-y
        const GInt nghbrmXmY                                 = nghbrmX != INVALID_CELLID ? tmpN(nghbrmX, 2) : -1;
        neighbor(cellId, cartesian::maxNoNghbrs<NDIM>() + 2) = nghbrmXmY;

        // -x+y
        const GInt nghbrmXpY                                 = nghbrmX != INVALID_CELLID ? tmpN(nghbrmX, 3) : -1;
        neighbor(cellId, cartesian::maxNoNghbrs<NDIM>() + 3) = nghbrmXpY;
      }

      // add tridiagonal nghbrs
      if(NDIM > 2) {
        TERMM(-1, "Not implemented!");
      }
    }
  }

  auto getCartesianGridData() const -> CartesianGridData<NDIM> { return CartesianGridData<NDIM>(*this); }

  auto totalSize() const -> GInt { return size() + m_noGhostsCells; }

 private:
  void setProperties() {
    for(GInt cellId = 0; cellId < noCells(); ++cellId) {
      const GBool isLeaf                     = noChildren(cellId) == 0;
      property(cellId, CellProperties::leaf) = isLeaf;
      m_noLeafCells += static_cast<GInt>(isLeaf);
    }
  };

  void determineBoundaryCells() {
    for(GInt cellId = 0; cellId < noCells(); ++cellId) {
      // is a partition cell determine for each if it can be a boundary cell (no existent parent)
      // parent has a cut with the boundary -> possible cut of child!
      if(parent(cellId) == -1 || property(parent(cellId), CellProperties::bndry)) {
        const GDouble cellLength                = lengthOnLvl(std::to_integer<GInt>(level(cellId)));
        property(cellId, CellProperties::bndry) = m_geometry->cutWithCell(center(cellId), cellLength);
        //        if(DEBUG_LEVEL > Debug_Level::min_debug && property(cellId, CellProperties::bndry)){
        if(property(cellId, CellProperties::bndry)) {
          GInt noNeighbors = 0;
          for(GInt nghbrId = 0; nghbrId < cartesian::maxNoNghbrs<NDIM>(); ++nghbrId) {
            if(neighbor(cellId, nghbrId) != INVALID_CELLID) {
              ++noNeighbors;
            }
          }
          if(cartesian::maxNoNghbrs<NDIM>() == noNeighbors) {
            //            cerr0 << "Removed boundary property cellId: " << cellId << " (" << center(cellId)[0] << ", " << center(cellId)[1]
            //                  << ") L:" << cellLength << std::endl;
            property(cellId, CellProperties::bndry) = false;
            logger << "Simplified bndry process!!!" << std::endl;
            //            TERMM(-1, "Cell marked as boundary, but is not on a boundary!");
          }
        }
      }
      m_noBndCells += static_cast<GInt>(property(cellId, CellProperties::bndry) && property(cellId, CellProperties::leaf));
    }
  }

#ifdef CLANG_COMPILER
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma ide diagnostic   ignored "cppcoreguidelines-pro-bounds-constant-array-index"
#endif
  void identifyBndrySurfaces() {
    if(m_axisAlignedBnd) {
      for(GInt surfId = 0; surfId < cartesian::maxNoNghbrs<NDIM>(); ++surfId) {
        m_bndrySurfaces.insert({static_cast<GString>(DirIdString[surfId]), Surface<DEBUG_LEVEL, NDIM>(this->getCartesianGridData())});
      }

      for(GInt cellId = 0; cellId < size(); ++cellId) {
        if(property(cellId, Cell::bndry)) {
          for(GInt dir = 0; dir < cartesian::maxNoNghbrs<NDIM>(); ++dir) {
            if(!hasNeighbor(cellId, dir)) {
              m_bndrySurfaces.at(static_cast<GString>(DirIdString[dir])).addCell(cellId, dir);
            }
          }
        }
      }
    } else {
      TERMM(-1, "Not implemented");
    }
  }
#ifdef CLANG_COMPILER
#pragma clang diagnostic pop
#endif

  void setupPeriodicConnections() {
    if(m_periodic) {
      std::vector<json>                    periodicBnds = get_all_items_with_value("periodic");
      std::unordered_map<GString, GString> periodicConnections;
      for(const auto& bnd : periodicBnds) {
        for(const auto& [surfName, props] : bnd.items()) {
          // check if periodic boundary should be handled as a boundary condition
          const auto generateBndry = config::opt_config_value<GBool>(props, "generateBndry", true);
          // skip if it should be
          if(!generateBndry) {
            periodicConnections.emplace(surfName, config::required_config_value<GString>(props, "connection"));
          }
        }
      }
      if(!periodicConnections.empty()) {
        logger << "Setting up periodic connections!" << std::endl;
        for(const auto& connection : periodicConnections) {
          if(periodicConnections.count(connection.second) != 0) {
            periodicConnections.erase(connection.second);
            addPeriodicConnection(bndrySurface(connection.first), bndrySurface(connection.second));
          }
        }
      }
    }
  }

  // todo: fix for refinement level changes
  // todo: simplify
  void addPeriodicConnection(const Surface<DEBUG_LEVEL, NDIM>& surfA, const Surface<DEBUG_LEVEL, NDIM>& surfB) {
    // connect cells of surfA and surfB
    for(const GInt cellIdA : surfA.getCellList()) {
      for(const GInt cellIdB : surfB.getCellList()) {
        GInt        notMatchingDir = -1;
        const auto& centerA        = center(cellIdA);
        const auto& centerB        = center(cellIdB);

        // find cells of surfA and surfB to connect
        for(GInt dir = 0; dir < NDIM; ++dir) {
          // connect cells that have one direction which is identical
          if(std::abs(centerA[dir] - centerB[dir]) > GDoubleEps) {
            if(notMatchingDir >= 0) {
              // periodic connection not possible since two directions don't match (3D)
              notMatchingDir = -1;
              break;
            }
            notMatchingDir = dir;
            continue;
          }
        }

        // cells need to be connected
        if(notMatchingDir >= 0) {
          // identify periodic direction for each cell
          const GInt nghbrDir = 2 * notMatchingDir;
          if(centerA[notMatchingDir] > centerB[notMatchingDir]) {
            if constexpr(DEBUG_LEVEL > Debug_Level::min_debug) {
              if(neighbor(cellIdB, nghbrDir) != INVALID_CELLID) {
                TERMM(-1, "Invalid set periodic connection! cellIdB:" + std::to_string(cellIdB) + " dir:" + std::to_string(nghbrDir));
              }
              if(neighbor(cellIdA, nghbrDir + 1) != INVALID_CELLID) {
                TERMM(-1, "Invalid set periodic connection! cellIdA:" + std::to_string(cellIdA) + " dir:" + std::to_string(nghbrDir + 1));
              }
            }
            neighbor(cellIdB, nghbrDir)     = cellIdA;
            neighbor(cellIdA, nghbrDir + 1) = cellIdB;
          } else {
            if constexpr(DEBUG_LEVEL > Debug_Level::min_debug) {
              if(neighbor(cellIdA, nghbrDir) != INVALID_CELLID) {
                TERMM(-1, "Invalid set periodic connection! cellIdA:" + std::to_string(cellIdA) + " dir:" + std::to_string(nghbrDir));
              }
              if(neighbor(cellIdB, nghbrDir + 1) != INVALID_CELLID) {
                TERMM(-1, "Invalid set periodic connection! cellIdB:" + std::to_string(cellIdB) + " dir:" + std::to_string(nghbrDir + 1));
              }
            }
            neighbor(cellIdA, nghbrDir)     = cellIdB;
            neighbor(cellIdB, nghbrDir + 1) = cellIdA;
          }
        }
      }
    }
  }


  void setWorkload() { TERMM(-1, "Not implemented!"); };
  void calculateOffspringsAndWeights() { TERMM(-1, "Not implemented!"); };

  void invalidate(const GInt begin, const GInt end) {
    std::fill(&parent(begin), &parent(end), INVALID_CELLID);
    std::fill(m_childIds.begin() + begin * cartesian::maxNoChildren<NDIM>(),
              m_childIds.begin() + end * cartesian::maxNoChildren<NDIM>(),
              INVALID_CELLID);
    std::fill(m_nghbrIds.begin() + begin * cartesian::maxNoNghbrsDiag<NDIM>(),
              m_nghbrIds.begin() + end * cartesian::maxNoNghbrsDiag<NDIM>(),
              INVALID_CELLID);
    std::fill(&globalId(begin), &globalId(end), INVALID_CELLID);
    std::fill(&level(begin), &level(end), std::byte(-1));
    std::fill(&center(begin), center(end), NAN);
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

  //  cartesian::Tree<DEBUG_LEVEL, NDIM> m_tree{};
  std::shared_ptr<GeometryManager<DEBUG_LEVEL, NDIM>> m_geometry;

  GInt m_noLeafCells   = 0;
  GInt m_noBndCells    = 0;
  GInt m_noGhostsCells = 0;

  GBool m_loadBalancing  = false;
  GBool m_diagonalNghbrs = false;
  GBool m_axisAlignedBnd = false;
  GBool m_periodic       = false;

  std::unordered_map<GString, Surface<DEBUG_LEVEL, NDIM>> m_bndrySurfaces;

  // Data containers
  std::vector<GInt> m_childIds{};
  std::vector<GInt> m_nghbrIds{};
  std::vector<GInt> m_noOffsprings{};

  std::vector<GFloat> m_weight{};
  std::vector<GFloat> m_workload{};
};

#endif // GRIDGENERATOR_CARTESIANGRID_H
