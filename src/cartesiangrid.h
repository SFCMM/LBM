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
#include "lbm/constants.h"
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
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::ref_currentHighestLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::partitionLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::ref_partitionLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::setMaxLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::boundingBox;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::setBoundingBox;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::maxLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::transformMaxLvl;

  CartesianGrid()                                        = default;
  ~CartesianGrid() override                              = default;
  CartesianGrid(const CartesianGrid&)                    = delete;
  CartesianGrid(CartesianGrid&&)                         = delete;
  auto operator=(const CartesianGrid&) -> CartesianGrid& = delete;
  auto operator=(CartesianGrid&&) -> CartesianGrid&      = delete;

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
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      if(m_bndrySurfaces.count(id) == 0) {
        TERMM(-1, "Invalid bndryId \"" + id + "\"");
      }
    }
    return m_bndrySurfaces.at(id);
  }

  auto hasBndrySurface(const GString& id) -> GBool { return m_bndrySurfaces.count(id) > 0; }


  /// Load the generated grid in-memory and set additional properties
  /// \param grid Generated grid.
  void loadGridInplace(const CartesianGridGen<DEBUG_LEVEL, NDIM>& grid, std::shared_ptr<ConfigurationAccess> properties) {
    m_config = properties;
    // grid.balance(); //todo: implement
    setCapacity(grid.capacity()); // todo: change for adaptation
    m_geometry              = grid.geometry();
    size()                  = grid.size();
    ref_currentHighestLvl() = grid.currentHighestLvl();
    ref_partitionLvl()      = grid.partitionLvl();
    setMaxLvl(grid.maxLvl());
    setBoundingBox(grid.boundingBox());
    transformMaxLvl(grid.lengthOnLvl(maxLvl()) / lengthOnLvl(maxLvl()));

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

    m_axisAlignedBnd = m_config->opt_config_value<GBool>("assumeAxisAligned", m_axisAlignedBnd);
    m_periodic       = m_config->has_any_key_value("type", "periodic");

    setProperties();

    determineBoundaryCells();
    identifyBndrySurfaces();
    setupPeriodicConnections();
    // todo: rename add interface cells...
    //    addGhostCells();
    addDiagonalNghbrs();

    //    for(auto& [name, srf] : m_bndrySurfaces) {
    //      srf.updateNeighbors();
    //    }
    //    if(m_loadBalancing) {
    //      setWorkload();
    //      calculateOffspringsAndWeights();
    //    }
  }

  /// Add ghost cells
  void addGhostCells() {
    struct PossibleBndGhost {
      std::vector<GInt>    dir;
      GBool                connected = false;
      GBool                invalid   = false;
      std::vector<GString> linkedSurfaces;
    };

    std::unordered_map<GInt, PossibleBndGhost> allPossibleBndryGhosts;
    std::set<GInt>                             boundaryCells;
    const GInt                                 bndryGhostOffset = size();

    // check all surfaces and add ghostcells in all missing dist directions
    for(const auto& [srfName, srf] : m_bndrySurfaces) {
      if(srf.hasBndryGhostCells()) {
        cerr0 << "srfName: " << srfName << " adds boundary ghost cells" << std::endl;
        logger << "srfName: " << srfName << " adds boundary ghost cells" << std::endl;
      }
      for(GInt cellId : srf.getCellList()) {
        // check all the main directions of the surface
        for(GInt nghbrDir = 0; nghbrDir < cartesian::maxNoNghbrs<NDIM>(); ++nghbrDir) {
          if(neighbor(cellId, nghbrDir) == INVALID_CELLID) {
            boundaryCells.emplace(cellId);

            if(srf.hasBndryGhostCells()) {
              // bndry ghost cell might need to be added
              allPossibleBndryGhosts[cellId].dir.emplace_back(nghbrDir);
              allPossibleBndryGhosts[cellId].linkedSurfaces.emplace_back(srfName);
              allPossibleBndryGhosts[cellId].connected = true;
            } else {
              allPossibleBndryGhosts[cellId].invalid = true;
            }
          }
        }
      }
    }

    // remove all invalid cells
    for(const auto& bndId : boundaryCells) {
      if(allPossibleBndryGhosts[bndId].invalid && !allPossibleBndryGhosts[bndId].connected) {
        allPossibleBndryGhosts.erase(bndId);
      }
    }

    cerr0 << "number of possible bndry ghost cells " << allPossibleBndryGhosts.size() << std::endl;

    auto addCell = [&](const GInt linkedCell, const GInt dir, const GString surfId) {
      const GInt ghostCellId                             = bndryGhostOffset + m_noGhostsCells;
      neighbor(linkedCell, dir)                          = ghostCellId;
      neighbor(ghostCellId, cartesian::oppositeDir(dir)) = linkedCell;

      GDouble length                               = 0.5 * lengthOnLvl(std::to_integer<GInt>(level(linkedCell)));
      level(ghostCellId)                           = std::byte(static_cast<GInt>(level(linkedCell)) + 1);
      center(ghostCellId)                          = center(linkedCell) + length * cartesian::dirVec<NDIM>(dir);
      property(ghostCellId, CellProperties::ghost) = true; // no valid solution
      property(ghostCellId, CellProperties::bndry) = true; // on boundary
      property(ghostCellId, CellProperties::solid) = true; // on solid side
      ++m_noGhostsCells;

      m_bndrySurfaces.at(surfId).addCell(ghostCellId, dir);
      return ghostCellId;
    };

    // generate cells
    for(auto& [linkedCell, value] : allPossibleBndryGhosts) {
      ASSERT(value.dir.size() <= 2, "Unsupported!" + std::to_string(value.dir.size()));

      const GString surfId = value.linkedSurfaces[0];
      //      GInt          ghostId = INVALID_CELLID;

      if(value.dir.size() == 1) {
        const GInt dir = value.dir[0];
        addCell(linkedCell, dir, surfId);
      } else {
        const VectorD<NDIM> surfNormalDir = m_bndrySurfaces.at(surfId).normal(linkedCell);

        for(GInt dir = 0; dir < cartesian::maxNoNghbrs<NDIM>(); ++dir) {
          const GDouble normalDir = surfNormalDir.dot(cartesian::dirVec<NDIM>(dir));
          if(normalDir > 0) {
            addCell(linkedCell, dir, surfId);
          }
        }

        //        const GInt diagonal = cartesian::inbetweenDiagDirs<NDIM>(value.dir[0], value.dir[1]);
        //        ghostId             = addCell(linkedCell, diagonal, surfId);
        //        value.dir.emplace_back(diagonal);

        for(GInt dir = 0; dir < cartesian::maxNoNghbrs<NDIM>(); ++dir) {
          if(neighbor(linkedCell, dir) == INVALID_CELLID) {
            const auto invalidDirId = std::find(value.dir.begin(), value.dir.end(), dir);
            if(invalidDirId != std::end(value.dir)) {
              value.dir.erase(invalidDirId);
            }
          }
        }
      }
    }

    // fix connections to neighbors on main directions
    for(const auto& [linkedCell, value] : allPossibleBndryGhosts) {
      const GString surfId = value.linkedSurfaces[0];
      for(GInt dirId = 0; dirId < value.dir.size(); ++dirId) {
        const GInt linkedDir = value.dir[dirId];
        const GInt ghostId   = neighbor(linkedCell, linkedDir);
        cerr0 << "ghost: " << ghostId << " is linked to " << linkedCell << " by " << linkedDir << std::endl;
        for(GInt dir = 0; dir < cartesian::maxNoNghbrs<NDIM>(); ++dir) {
          if(neighbor(ghostId, dir) == INVALID_CELLID) {
            const GInt linkedNghbrId = neighbor(linkedCell, dir);
            const GInt nghbrId       = (linkedNghbrId != INVALID_CELLID) ? neighbor(linkedNghbrId, linkedDir) : INVALID_CELLID;
            if(nghbrId != INVALID_CELLID && nghbrId != ghostId) {
              neighbor(ghostId, dir)                         = nghbrId;
              neighbor(nghbrId, cartesian::oppositeDir(dir)) = ghostId;
            }
          }
        }
      }
      m_bndrySurfaces.at(surfId).removeCell(linkedCell);
    }

    //    addDiagonalNghbrs(bndryGhostOffset);

    if(m_noGhostsCells > 0) {
      cerr0 << "Added ghostlayers no cells:" << m_noGhostsCells << std::endl;
      logger << "Added ghostlayers no cells:" << m_noGhostsCells << std::endl;
    }
  }

  /// Add the diagonal(2D/3D) and/or tridiagonal (3D) to the neighbor connections of each cell.
  void addDiagonalNghbrs(const GInt offset = 0) {
    m_diagonalNghbrs = true;
    auto tmpNghbr    = m_nghbrIds;

    auto tmpN = [&](const GInt cellId, const GInt dir) { return tmpNghbr[cellId * cartesian::maxNoNghbrs<NDIM>() + dir]; };


    for(GInt cellId = offset; cellId < size() + m_noGhostsCells; ++cellId) {
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

  auto getCartesianGridData() -> CartesianGridData<NDIM> { return CartesianGridData<NDIM>(*this); }

  auto totalSize() const -> GInt { return size() + m_noGhostsCells; }

  auto isBndryCell(const GInt cellId) const -> GBool { return property(cellId, CellProperties::bndry); }

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
      // is a partition cell determine for each if it can be a boundary cell
      // cell has no parent -> might have cut
      // parent has a cut with the boundary -> possible cut of child!
      if(parent(cellId) == -1 || property(parent(cellId), CellProperties::bndry)) {
        const GDouble cellLength = lengthOnLvl(std::to_integer<GInt>(level(cellId)));

        // check for cut with geometry
        property(cellId, CellProperties::bndry) = m_geometry->cutWithCell(center(cellId), cellLength);
        //        if(DEBUG_LEVEL > Debug_Level::min_debug && property(cellId, CellProperties::bndry)){

        // we currently only care for cells which have missing neighbors!
        if(property(cellId, CellProperties::bndry)) {
          GInt noNeighbors = 0;
          for(GInt nghbrId = 0; nghbrId < cartesian::maxNoNghbrs<NDIM>(); ++nghbrId) {
            if(neighbor(cellId, nghbrId) != INVALID_CELLID) {
              ++noNeighbors;
            }
          }

          // remove cells that have no missing neighbor
          if(cartesian::maxNoNghbrs<NDIM>() == noNeighbors) {
            //            cerr0 << "Removed boundary property cellId: " << cellId << " (" << center(cellId)[0] << ", " <<
            //            center(cellId)[1]
            //                  << ") L:" << cellLength << std::endl;
            property(cellId, CellProperties::bndry) = false;
            if(property(cellId, CellProperties::leaf)) {
              logger << "Simplified bndry process!!! Cell had cut with surface but is not missing Neighbors." << std::endl;
            }
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
    json bndryConfig = m_config->getObject("boundary");

    cerr0 << "Create boundary surfaces for " << m_geometry->noObjects() << " geometries" << std::endl;

    std::set<std::pair<GInt, GInt>> assignedBnds;

    // iterate over all geometries
    for(const auto& [surfName, surfConfig] : bndryConfig.items()) {
      const GInt noBnds = surfConfig.size();
      for(const auto& [surfDirName, config] : surfConfig.items()) {
        const GString surfNameAp = (noBnds > 1) ? surfName + "_" + surfDirName : surfName;
        cerr0 << "Surface name: " << surfNameAp << std::endl;
        cerr0 << "bndConfig: " << config << std::endl;
        m_bndrySurfaces.insert(std::make_pair(surfNameAp, Surface<DEBUG_LEVEL, NDIM>(this->getCartesianGridData(), &property(0))));

        // use "all" to set all direction for this bnd
        const GInt dirBegin = surfDirName == "all" ? 0 : dirIdString2Id(surfDirName);
        const GInt dirEnd   = surfDirName == "all" ? cartesian::maxNoNghbrs<NDIM>() : dirIdString2Id(surfDirName) + 1;

        for(GInt dir = dirBegin; dir < dirEnd; ++dir) {
          for(GInt cellId = 0; cellId < size(); ++cellId) {
            if(property(cellId, Cell::bndry)) {
              const GDouble cellLength = lengthOnLvl(std::to_integer<GInt>(level(cellId)));
              if(!hasNeighbor(cellId, dir)) {
                // cell has cut with the boundary surface
                if(m_geometry->cutWithCell(surfName, center(cellId), cellLength)) {
                  const auto [iter, added] = assignedBnds.insert({cellId, dir});

                  if(added) {
                    // todo: diagonal missing cells are not assigned!
                    //                    cerr0 << surfName << " : " << cellId << std::endl;
                    m_bndrySurfaces.at(surfNameAp).addCell(cellId, dir);
                  } else {
                    cerr0 << "cellId: " << cellId << " was already assigned a bnd in direction " << dir << std::endl;
                  }
                }
              }
            }
          }
        }
        if(m_bndrySurfaces.at(surfNameAp).size() == 0) {
          //            m_bndrySurfaces.erase(surfNameAp);
          cerr0 << "WARNING: surface " << surfNameAp << " has no cells!" << std::endl;
          logger << "WARNING: surface " << surfNameAp << " has no cells!" << std::endl;
        } else {
          cerr0 << "Surface assigned " << m_bndrySurfaces.at(surfNameAp).size() << " cells" << std::endl;
        }
      }
    }
  }
#ifdef CLANG_COMPILER
#pragma clang diagnostic pop
#endif

  void setupPeriodicConnections() {
    if(m_periodic) {
      json                                 bndryConfig = m_config->getObject("boundary");
      std::unordered_map<GString, GString> periodicConnections;

      for(const auto& [geometryName, geomConf] : bndryConfig.items()) {
        std::vector<json> periodicBnds = m_config->get_all_items_with_value("periodic");

        for(const auto& [surfName, surfConf] : geomConf.items()) {
          if(config::opt_config_value(surfConf, "type", static_cast<GString>("notset")) == "periodic") {
            // check if periodic boundary should be handled as a boundary condition
            const auto generateBndry = config::opt_config_value<GBool>(surfConf, "generateBndry", true);
            // skip if it should be
            if(!generateBndry) {
              periodicConnections.emplace(geometryName + "_" + surfName, config::required_config_value<GString>(surfConf, "connection"));
              cerr0 << "Add connection: " << geometryName + "_" + surfName << " with " << surfConf << std::endl;
            }
          }
        }
      }

      if(!periodicConnections.empty()) {
        logger << "Setting up periodic connections!" << std::endl;
        for(const auto& connection : periodicConnections) {
          if(periodicConnections.count(connection.second) != 0) {
            periodicConnections.erase(connection.second);
            addPeriodicConnection(bndrySurface(connection.first), bndrySurface(connection.second));
          } else {
            TERMM(-1, "Invalid periodic setup!");
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
            if constexpr(DEBUG_LEVEL >= Debug_Level::debug) {
              if(neighbor(cellIdB, nghbrDir) != INVALID_CELLID) {
                TERMM(-1, "Invalid set periodic connection! cellIdB:" + std::to_string(cellIdB) + " dir:" + std::to_string(nghbrDir));
              }
              if(neighbor(cellIdA, nghbrDir + 1) != INVALID_CELLID) {
                TERMM(-1, "Invalid set periodic connection! cellIdA:" + std::to_string(cellIdA) + " dir:" + std::to_string(nghbrDir + 1));
              }
            }
            neighbor(cellIdB, nghbrDir)     = cellIdA;
            neighbor(cellIdA, nghbrDir + 1) = cellIdB;
            if constexpr(DEBUG_LEVEL == Debug_Level::max_debug) {
              logger << "connected " << cellIdA << " with " << cellIdB << std::endl;
            }
          } else {
            if constexpr(DEBUG_LEVEL >= Debug_Level::debug) {
              if(neighbor(cellIdA, nghbrDir) != INVALID_CELLID) {
                TERMM(-1, "Invalid set periodic connection! cellIdA:" + std::to_string(cellIdA) + " dir:" + std::to_string(nghbrDir));
              }
              if(neighbor(cellIdB, nghbrDir + 1) != INVALID_CELLID) {
                TERMM(-1, "Invalid set periodic connection! cellIdB:" + std::to_string(cellIdB) + " dir:" + std::to_string(nghbrDir + 1));
              }
            }
            neighbor(cellIdA, nghbrDir)     = cellIdB;
            neighbor(cellIdB, nghbrDir + 1) = cellIdA;

            if constexpr(DEBUG_LEVEL == Debug_Level::max_debug) {
              logger << "connected " << cellIdA << " with " << cellIdB << std::endl;
            }
          }
        } else {
          logger << "No periodic connection found for cells: " << cellIdA << " and " << cellIdB << std::endl;
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

  std::shared_ptr<ConfigurationAccess> m_config;
};

#endif // GRIDGENERATOR_CARTESIANGRID_H
