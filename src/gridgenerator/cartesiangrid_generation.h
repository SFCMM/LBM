#ifndef GRIDGENERATOR_CARTESIANGRID_GENERATION_H
#define GRIDGENERATOR_CARTESIANGRID_GENERATION_H
#include <sfcmm_common.h>
#include "cartesiangrid_base.h"

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class CartesianGridGen : public BaseCartesianGrid<DEBUG_LEVEL, NDIM> {
 public:
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::partitionLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::maxLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::lengthOnLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::cog;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::increaseCurrentHighestLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::currentHighestLvl;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::geometry;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::property;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::parent;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::level;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::globalId;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::center;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::empty;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::size;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::capacity;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::boundingBox;
  using BaseCartesianGrid<DEBUG_LEVEL, NDIM>::transformMaxLvl;

  using PropertyBitsetType = grid::cell::BitsetType;
  using ChildListType      = std::array<GInt, cartesian::maxNoChildren<NDIM>()>;

  explicit CartesianGridGen(const GInt maxNoCells) { setCapacity(maxNoCells); }
  ~CartesianGridGen() override              = default;
  CartesianGridGen(const CartesianGridGen&) = delete;
  CartesianGridGen(CartesianGridGen&&)      = delete;
  auto operator=(const CartesianGridGen&) -> CartesianGridGen& = delete;
  auto operator=(CartesianGridGen&&) -> CartesianGridGen& = delete;

  void setCapacity(const GInt _capacity) override {
    if(!empty()) {
      TERMM(-1, "Invalid operation tree already allocated.");
    }
    if(_capacity < capacity()) {
      logger << "WARNING: memory allocation requested. But capacity is sufficient CartesianGridGen::setCapacity() requested " << capacity()
             << " but allocated " << capacity() << std::endl;
      return;
    }
    m_noChildren.resize(_capacity);
    m_nghbrIds.resize(_capacity);
    m_childIds.resize(_capacity);
    m_rfnDistance.resize(_capacity);
    BaseCartesianGrid<DEBUG_LEVEL, NDIM>::setCapacity(_capacity);
  }

  void reset() override {
    m_noChildren.clear();
    m_nghbrIds.clear();
    m_childIds.clear();
    m_rfnDistance.clear();
    BaseCartesianGrid<DEBUG_LEVEL, NDIM>::clear();
  }

  void setMaxLvl(const GInt _maxLvl) override {
    m_levelOffsets.resize(_maxLvl + 1);
    BaseCartesianGrid<DEBUG_LEVEL, NDIM>::setMaxLvl(_maxLvl);
  }

  /// Create the grid that is used for partitioning. This grid has the level of the option provided in the grid
  /// configuration file. The grid up to this level is always produced on a single MPI rank.
  /// \param partitioningLvl Level of the partitioning grid.
  void createPartitioningGrid(const GInt partitioningLvl) {
    RECORD_TIMER_START(TimeKeeper[Timers::GridPart]);
    if(capacity() < 1) {
      TERMM(-1, "Invalid grid capacity.");
    }
    if(geometry() == nullptr) {
      TERMM(-1, "No geometry set");
    }

    partitionLvl() = partitioningLvl;
    logger << SP1 << "Create partitioning grid with level " << partitionLvl() << std::endl;
    std::cout << SP1 << "Create partitioning grid with level " << partitionLvl() << std::endl;

    logger << SP2 << "+ initial cube length: " << lengthOnLvl(0) << std::endl;
    std::cout << SP2 << "+ initial cube length: " << lengthOnLvl(0) << std::endl;

    // make sure we have set some level...
    if(partitioningLvl > maxLvl()) {
      cerr0 << "WARNING: No maximum level set -> set to " << partitioningLvl << std::endl;
      logger << "WARNING: No maximum level set -> set to " << partitioningLvl << std::endl;
      setMaxLvl(partitioningLvl);
    }

    // use lazy initialization for grid generation and make sure final partitioning grid starts in the beginning
    if(isEven(partitionLvl())) {
      // initial cell placed in the beginning
      m_levelOffsets[0] = {0, 1};
      m_levelOffsets[1] = {capacity() - cartesian::maxNoChildren<NDIM>(), capacity()};
    } else {
      // initial cell placed at the end
      m_levelOffsets[0] = {capacity() - 1, capacity()};
      m_levelOffsets[1] = {0, cartesian::maxNoChildren<NDIM>()};
    }

    const GInt begin                       = m_levelOffsets[0].begin;
    center(begin)                          = Point<NDIM>(cog().data());
    globalId(begin)                        = begin;
    property(begin, CellProperties::bndry) = true;
    size()                                 = 1;

    //  Refine to min level
    for(GInt l = 0; l < partitionLvl(); l++) {
      const GInt prevLevelBegin   = m_levelOffsets[l].begin;
      const GInt prevLevelEnd     = m_levelOffsets[l].end;
      const GInt prevLevelNoCells = prevLevelEnd - prevLevelBegin;
      if(DEBUG_LEVEL > Debug_Level::no_debug) {
        logger << "LevelOffset " << l << ":[" << prevLevelBegin << ", " << prevLevelEnd << "]" << std::endl;
      }


      if(m_levelOffsets[l].begin == 0) {
        // m_capacity - (prevLevelNoCells) * cartesian::maxNoChildren<NDIM>()
        // from the end - maximum number of cells at the current level
        const GInt newLevelBegin = capacity() - (prevLevelNoCells)*cartesian::maxNoChildren<NDIM>();
        m_levelOffsets[l + 1]    = {newLevelBegin, capacity()};

        if(prevLevelEnd > newLevelBegin) {
          outOfMemory(l + 1);
        }
      } else {
        const GInt newLevelEnd = (prevLevelNoCells)*cartesian::maxNoChildren<NDIM>();
        // from the start to the maximum number of cells at the current level
        m_levelOffsets.at(l + 1) = {0, newLevelEnd};

        if(prevLevelBegin < newLevelEnd) {
          outOfMemory(l + 1);
        }
      }

      refineGrid(l);
    }

    std::fill(&parent(0), &parent(capacity() - 1), INVALID_CELLID);
    reorderHilbertCurve();

    RECORD_TIMER_STOP(TimeKeeper[Timers::GridPart]);
  }

  /// Uniformly refine the grid up to the provided level.
  /// \param uniformLvl Level of uniform refinement.
  void uniformRefineGrid(const GInt uniformLevel) {
    RECORD_TIMER_START(TimeKeeper[Timers::GridUniform]);
    logger << SP1 << "Uniformly refine grid to level " << uniformLevel << std::endl;
    std::cout << SP1 << "Uniformly refine grid to level " << uniformLevel << std::endl;

    if(partitionLvl() == uniformLevel) {
      return;
    }

    for(GInt lvl = partitionLvl(); lvl < uniformLevel; lvl++) {
      m_levelOffsets[lvl + 1] = {size(), size() + levelSize(m_levelOffsets[lvl]) * cartesian::maxNoChildren<NDIM>()};
      if(m_levelOffsets[lvl + 1].end > capacity()) {
        outOfMemory(lvl + 1);
      }

      refineGrid<true>(lvl);
    }
    RECORD_TIMER_STOP(TimeKeeper[Timers::GridUniform]);
  }

  /// Refine the cells that have been marked for refinement.
  /// \param noCellsToRefine The number of cells that have been marked.
  void refineMarkedCells(const GInt noCellsToRefine) {
    if(noCellsToRefine == 0) {
      logger << "WARNING: refineMarkedCells called but nothing to do" << std::endl;
      return;
    }

    logger << SP1 << "Refining marked cells to level " << currentHighestLvl() + 1 << std::endl;
    std::cout << SP1 << "Refining marked cells to level " << currentHighestLvl() + 1 << std::endl;
    // update the offsets
    m_levelOffsets[currentHighestLvl() + 1] = {size(), size() + noCellsToRefine * cartesian::maxNoChildren<NDIM>()};
    if(m_levelOffsets[currentHighestLvl() + 1].end > capacity()) {
      outOfMemory(currentHighestLvl() + 1);
    }
    logger << SP2 << "* cells to refine: " << noCellsToRefine << std::endl;
    std::cout << SP2 << "* cells to refine: " << noCellsToRefine << std::endl;


    refineGrid<true, false>(currentHighestLvl());
  }

  /// Mark all boundary cells for refinement
  /// \return Number of cells marked for refinement
  auto markBndryCells() -> GInt {
    logger << SP2 << "* marking bndry cells " << std::endl;
    std::cout << SP2 << "* marking bndry cells " << std::endl;
    GInt markedCells = 0;
    for(GInt cellId = m_levelOffsets[currentHighestLvl()].begin; cellId < m_levelOffsets[currentHighestLvl()].end; ++cellId) {
      if(property(cellId, CellProperties::bndry)) {
        property(cellId, CellProperties::toRefine) = true;
        markedCells++;
      }
    }
    return markedCells;
  }

  void save(const GString& fileName, const json& gridOutConfig) const override {
    // Grid output configuration
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    GString              filter    = config::opt_config_value(gridOutConfig, "cellFilter", GString("leafCells"));
    GString              format    = config::opt_config_value(gridOutConfig, "format", GString("ASCII"));
    GString              type      = config::opt_config_value(gridOutConfig, "type", GString("point"));
    std::vector<GString> outvalues = config::opt_config_value(gridOutConfig, "outputValues", std::vector<GString>({"level"}));
    GInt                 outputLvl = config::opt_config_value(gridOutConfig, "outputLvl", 0);

    // Cell filter functions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // only output the lowest level
    std::function<GBool(GInt)> isLowestLevel = [&](GInt cellId) { return std::to_integer<GInt>(level(cellId)) == partitionLvl(); };

    // only output the highest level
    std::function<GBool(GInt)> isHighestLevel = [&](GInt cellId) { return std::to_integer<GInt>(level(cellId)) == currentHighestLvl(); };

    // only output leaf cells (i.e. cells without children)
    std::function<GBool(GInt)> isLeaf = [&](GInt cellId) { return m_noChildren[cellId] == 0; };

    // only output the given level
    std::function<GBool(GInt)> isTargetLevel = [&](GInt cellId) { return std::to_integer<GInt>(level(cellId)) == outputLvl; };


    std::function<GBool(GInt)>& outputFilter = isLeaf;
    if(filter == "highestLvl") {
      outputFilter = isHighestLevel;
    } else if(filter == "lowestLvl" || filter == "partitionLvl") {
      outputFilter = isLowestLevel;
    } else if(filter == "leafCells") {
      outputFilter = isLeaf;
    } else if(filter == "targetLvl") {
      if(!config::has_config_value(gridOutConfig, "outputLvl")) {
        TERMM(-1, "Required value not set!");
      }
      outputFilter = isTargetLevel;
      if(outputLvl < partitionLvl()) {
        TERMM(-1, "Outputting a lvl below the partition lvl is not possible!");
      }
    } else {
      TERMM(-1, "Unknown output filter " + filter);
    }

    std::vector<IOIndex>              index;
    std::vector<std::vector<GString>> values;
    cerr0 << "Selected output values:";
    for(const auto& outputvalue : outvalues) {
      if(outputvalue == "level") {
        // todo:replace type
        index.emplace_back(IOIndex{"Level", "int64"});
        values.emplace_back(toStringVector(level(), size()));
        cerr0 << " level ";
      } else if(outputvalue == "noChildren") {
        // todo:replace type
        index.emplace_back(IOIndex{"NoChildren", "int64"});
        values.emplace_back(toStringVector(m_noChildren, size()));
        cerr0 << " noChildren ";
      } else {
        logger << "WARNING: The output value " + outputvalue + " is not a valid output!" << std::endl;
      }
    }
    cerr0 << std::endl;

    if(format == "ASCII") {
      ASCII::writePointsCSV<NDIM>(fileName, size(), center(), index, values, outputFilter);
    } else if(format == "VTK") {
      VTK::ASCII::writePoints<NDIM>(fileName, size(), center(), index, values, outputFilter);
    } else if(format == "VTKB") {
      // todo: rename format
      VTK::BINARY::writePoints<NDIM>(fileName, size(), center(), index, values, outputFilter);
    } else {
      TERMM(-1, "Unknown output format " + format);
    }
  }

  /// Move the leaf nodes to the surface of the geometry.
  void transformMaxRfnmtLvlToExtent() {
    BoundingBoxCT<NDIM> actualExtent;
    for(GInt dir = 0; dir < NDIM; ++dir) {
      actualExtent.min(dir) = std::numeric_limits<GDouble>::max();
      actualExtent.max(dir) = std::numeric_limits<GDouble>::min();
    }

    for(GInt cellId = m_levelOffsets[maxLvl()].begin; cellId < m_levelOffsets[maxLvl()].end; ++cellId) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        if(actualExtent.min(dir) > center(cellId, dir)) {
          actualExtent.min(dir) = center(cellId, dir);
        }
        if(actualExtent.max(dir) < center(cellId, dir)) {
          actualExtent.max(dir) = center(cellId, dir);
        }
      }
    }

    cerr0 << "actual nodal extent: " << actualExtent.str() << std::endl;
    GInt    alignDir            = 1;
    GDouble transformationValue = boundingBox().max(alignDir) / (actualExtent.max(alignDir) - actualExtent.min(alignDir));


    for(GInt cellId = m_levelOffsets[maxLvl()].begin; cellId < m_levelOffsets[maxLvl()].end; ++cellId) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        if(dir == alignDir) {
          center(cellId, dir) =
              center(cellId, dir) * transformationValue - (transformationValue * actualExtent.max(dir) - boundingBox().max(dir));
        } else {
          center(cellId, dir) = center(cellId, dir) * transformationValue;
        }
      }
    }
    transformMaxLvl(transformationValue);

    // afterwards there are cells which will be outside
    deleteOutsideCells<true>(maxLvl());
  }

  static constexpr auto memorySizePerCell() -> GInt {
    return sizeof(GInt) * (1 + 1 + 1 + 1 + 2) // m_parentId, m_globalId, m_noChildren, m_rfnDistance, m_levelOffsets
           + sizeof(Point<NDIM>)              // m_center
           + sizeof(NeighborList<NDIM>)       // m_nghbrIds
           + sizeof(ChildList<NDIM>)          // m_childIds
           + sizeof(PropertyBitsetType)       // m_properties
           + 1;                               // m_level
  }

  [[nodiscard]] auto child(const GInt id, const GInt childId) const -> GInt { return m_childIds[id].c[childId]; }

  [[nodiscard]] auto neighbor(const GInt id, const GInt dir) const -> GInt override { return m_nghbrIds[id].n[dir]; }

 protected:
  [[nodiscard]] auto neighbor(const GInt id, const GInt dir) -> GInt& { return m_nghbrIds[id].n[dir]; }


 private:
  void outOfMemory(GInt _level) {
    cerr0 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    logger << "ERROR: Not enough memory to generate grid! Increase maxNoCells: " << capacity() << std::endl;
    cerr0 << "ERROR: Not enough memory to generate grid! Increase maxNoCells: " << capacity() << std::endl;
    logger << "level " << _level - 1 << " [" << m_levelOffsets[_level - 1].begin << ", " << m_levelOffsets[_level - 1].end << "]"
           << std::endl;
    cerr0 << "level " << _level - 1 << " [" << m_levelOffsets[_level - 1].begin << ", " << m_levelOffsets[_level - 1].end << "]"
          << std::endl;
    logger << "level " << _level << " [" << m_levelOffsets[_level].begin << ", " << m_levelOffsets[_level].end << "]" << std::endl;
    cerr0 << "level " << _level << " [" << m_levelOffsets[_level].begin << ", " << m_levelOffsets[_level].end << "]" << std::endl;
    cerr0 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

    TERMM(-1, "Out of memory!");
  }

  /// Refine grid to an one higher level.
  /// \tparam DISTRIBUTED Grid is using MPI
  /// \tparam UNIFORM Grid is uniform
  /// \param lvlToBeRefined Level to be refined
  template <GBool DISTRIBUTED = false, GBool UNIFORM = true>
  void refineGrid(const GInt lvlToBeRefined) {
    if(DISTRIBUTED && !MPI::isSerial()) {
      // updateHaloOffsets();
    }

    if(UNIFORM) {
      refineGrid(m_levelOffsets, lvlToBeRefined);
      if(DISTRIBUTED && !MPI::isSerial()) {
        // todo:implement
        // refineGrid(m_haloOffsets, l);
      }
    } else {
      refineGridMarkedOnly(m_levelOffsets, lvlToBeRefined);
      if(DISTRIBUTED && !MPI::isSerial()) {
        // todo:implement
        // refineGridMarkedOnly(m_haloOffsets, l);
      }
    }
    size() = m_levelOffsets[lvlToBeRefined + 1].end;

    findChildLevelNghbrs(m_levelOffsets, lvlToBeRefined);
    if(DISTRIBUTED && !MPI::isSerial()) {
      // todo:implement
      // findChildLevelNeighbors(m_haloOffsets, l);
    }

    if(DISTRIBUTED && !MPI::isSerial()) {
      // todo:implement
      //        deleteOutsideCellsParallel(l + 1);
    } else {
      deleteOutsideCells(lvlToBeRefined + 1);
    }

    size() = m_levelOffsets[lvlToBeRefined + 1].end;
    increaseCurrentHighestLvl();
    logger.updateAttributes();
  }

  void refineGrid(const std::vector<LevelOffsetType>& levelOffset, const GInt _level) {
    logger << SP2 << "+ refining grid on level: " << _level << std::endl;
    std::cout << SP2 << "+ refining grid on level: " << _level << std::endl;

    ASSERT(_level + 1 <= maxLvl(), "Invalid refinement level! " + std::to_string(_level + 1) + ">" + std::to_string(maxLvl()));

    // refine all cells on the given level
#ifdef _OPENMP
#pragma omp parallel default(none) shared(levelOffset, _level)
    {
#endif

      const GInt firstCellOfLvl = levelOffset[_level].begin;
      const GInt lastCellOfLvl  = levelOffset[_level].end;
      const GInt refinedLvl     = _level + 1;
#ifdef _OPENMP
#pragma omp for
#endif
      for(GInt cellId = firstCellOfLvl; cellId < lastCellOfLvl; ++cellId) {
        const GInt cellCount = cellId - firstCellOfLvl;
        refineCell(cellId, levelOffset[refinedLvl].begin + cellCount * cartesian::maxNoChildren<NDIM>());
      }
#ifdef _OPENMP
    }
#endif
  }

  void refineGridMarkedOnly(const std::vector<LevelOffsetType>& levelOffset, const GInt _level) {
    GInt refinedCells = 0;
    for(GInt cellId = levelOffset[_level].begin; cellId < levelOffset[_level].end; ++cellId) {
      if(property(cellId, CellProperties::toRefine)) {
        refineCell(cellId, levelOffset[_level + 1].begin + refinedCells * cartesian::maxNoChildren<NDIM>());
        ++refinedCells;
      }
    }
  }

  void refineCell(const GInt cellId, const GInt offset) {
    if(DEBUG_LEVEL > Debug_Level::debug) {
      logger << "refine cell " << cellId << " with offset " << offset << std::endl;
    }
    const GInt    refinedLvl       = std::to_integer<GInt>(level(cellId)) + 1;
    const GDouble refinedLvlLength = lengthOnLvl(refinedLvl);

    for(GInt childId = 0; childId < cartesian::maxNoChildren<NDIM>(); ++childId) {
      const GInt childCellId = offset + childId;
      center(childCellId) =
          center(cellId)
          + HALF * refinedLvlLength
                * Point<NDIM>(cartesian::childDir[childId].data()); // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
      level(childCellId)    = static_cast<std::byte>(refinedLvl);
      parent(childCellId)   = cellId;
      globalId(childCellId) = childCellId;

      // reset since we overwrite previous levels
      m_noChildren[childCellId] = 0;
      property(childCellId).reset();
      m_childIds[childCellId] = {INVALID_LIST<cartesian::maxNoChildren<NDIM>()>()};
      m_nghbrIds[childCellId] = {INVALID_LIST<cartesian::maxNoNghbrs<NDIM>()>()};

      // if parent is a boundary cell check for children as well
      if(property(cellId, CellProperties::bndry)) {
        property(childCellId, CellProperties::bndry) = cellHasCut(childCellId, refinedLvlLength);
      }

      // update parent
      m_childIds[cellId].c[childId] = childCellId;
      m_noChildren[cellId]++;
    }
  }

  void findChildLevelNghbrs(const std::vector<LevelOffsetType>& levelOffset, const GInt _level) {
    // check all children at the given level
    for(GInt parentId = levelOffset[_level].begin; parentId < levelOffset[_level].end; ++parentId) {
      const GInt* __restrict children = &m_childIds[parentId].c[0];
      for(GInt childId = 0; childId < cartesian::maxNoChildren<NDIM>(); ++childId) {
        if(children[childId] == INVALID_CELLID) {
          // no child
          continue;
        }

        // check all neighbors
        GInt* __restrict neighbors = &m_nghbrIds[children[childId]].n[0];
        for(GInt dir = 0; dir < cartesian::maxNoNghbrs<NDIM>(); ++dir) {
          // neighbor direction not set
          if(neighbors[dir] == INVALID_CELLID) {
            const GInt nghbrId =
                static_cast<GInt>(cartesian::nghbrInside[childId][dir]); // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
            // neighbor is within the same parent cell
            if(nghbrId != INVALID_CELLID) {
              neighbors[dir] = children[nghbrId];
            } else {
              const GInt parentLvlNeighborChildId = static_cast<GInt>(
                  cartesian::nghbrParentChildId[childId][dir]); // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
              ASSERT(parentLvlNeighborChildId > INVALID_CELLID, "The definition of nghbrParentChildId is wrong! "
                                                                "childId: "
                                                                    + std::to_string(childId) + " dir " + std::to_string(dir));

              const GInt parentLvlNghbrId = m_nghbrIds[parentId].n[dir];
              if(parentLvlNghbrId != INVALID_CELLID && parentLvlNeighborChildId != INVALID_CELLID
                 && m_childIds[parentLvlNghbrId].c[parentLvlNeighborChildId] != INVALID_CELLID) {
                neighbors[dir] = m_childIds[parentLvlNghbrId].c[parentLvlNeighborChildId];
              }
            }

            if(DEBUG_LEVEL >= Debug_Level::debug && neighbors[dir] != INVALID_CELLID
               && (center(neighbors[dir]) - center(children[childId])).norm()
                      > 1.9 * lengthOnLvl(std::to_integer<GInt>(level(children[childId])))) {
              cerr0 << "neighbors[dir] " << neighbors[dir] << " cellId " << children[childId] << std::endl;
              cerr0 << "neighbors " << strStreamify<NDIM>(center(neighbors[dir])).str() << std::endl;
              cerr0 << "neighbors " << strStreamify<NDIM>(center(children[childId])).str() << std::endl;
              cerr0 << "ndiff " << (center(neighbors[dir]) - center(children[childId])).norm() << " vs "
                    << lengthOnLvl(std::to_integer<GInt>(level(children[childId]))) << std::endl;
              cerr0 << "parentId " << parentId << " np " << m_nghbrIds[parentId].n[dir] << std::endl;
              cerr0 << "parent " << strStreamify<NDIM>(center(parentId)).str() << std::endl;
              cerr0 << "pdiff " << (center(parentId) - center(m_nghbrIds[parentId].n[dir])).norm() << " vs "
                    << lengthOnLvl(std::to_integer<GInt>(level(parentId))) << std::endl;

              TERMM(-1, "Invalid neighbor");
            }
          }
        }
      }
    }
  }

  template <GBool CHECKALL = false>
  void deleteOutsideCells(const GInt _level) {
    markOutsideCells<CHECKALL>(m_levelOffsets, _level);

    // delete cells that have been marked as being outside
    for(GInt cellId = m_levelOffsets[_level].end - 1; cellId >= m_levelOffsets[_level].begin; --cellId) {
      ASSERT(!property(cellId, CellProperties::bndry)
                 || property(cellId, CellProperties::inside) == property(cellId, CellProperties::bndry),
             "Properties not set correctly! bndry implies IsInside!");
      // remove cell since it is not inside
      if(!property(cellId, CellProperties::inside)) {
        const GInt parentId = parent(cellId);

        // partitionlvl doesn't have parents
        if(_level != partitionLvl()) {
          ASSERT(parentId != INVALID_CELLID, "Invalid parentId! (cellId: " + std::to_string(cellId) + ")");

          // remove from parent
          updateParent(parent(cellId), cellId, INVALID_CELLID);
          --m_noChildren[parentId];
        }

        // remove from neighbors
        for(GInt dir = 0; dir < cartesian::maxNoNghbrs<NDIM>(); ++dir) {
          const GInt nghbrCellId = m_nghbrIds[cellId].n[dir];
          if(nghbrCellId != INVALID_CELLID) {
            m_nghbrIds[nghbrCellId].n[cartesian::oppositeDir(dir)] = INVALID_CELLID;
          }
        }
        if(cellId != m_levelOffsets[_level].end - 1) {
          // copy an inside cell to the current position to fill the hole
          copyCell(m_levelOffsets[_level].end - 1, cellId);
        }
        m_levelOffsets[_level].end--;
      }
    }
    size() = levelSize(m_levelOffsets[_level]);
    logger << SP3 << "* grid has " << size() << " cells" << std::endl;
    std::cout << SP3 << "* grid has " << size() << " cells" << std::endl;
  }

  template <GBool CHECKALL = false>
  void markOutsideCells(const std::vector<LevelOffsetType>& levelOffset, const GInt _level) {
    if(CHECKALL) {
      for(GInt cellId = levelOffset[_level].begin; cellId < levelOffset[_level].end; ++cellId) {
        const GBool isBndryCell                  = cellHasCut(cellId);
        property(cellId, CellProperties::bndry)  = isBndryCell;
        property(cellId, CellProperties::inside) = isBndryCell || pointIsInside(center(cellId));
      }
    } else {
      // reset marked property
      for(GInt cellId = levelOffset[_level].begin; cellId < levelOffset[_level].end; ++cellId) {
        property(cellId, CellProperties::marked) = false;
      }

      for(GInt cellId = levelOffset[_level].begin; cellId < levelOffset[_level].end; ++cellId) {
        if(property(cellId, CellProperties::marked)) {
          continue;
        }
        property(cellId, CellProperties::marked) = true;
        const GBool isBndryCell                  = property(cellId, CellProperties::bndry);
        property(cellId, CellProperties::inside) = isBndryCell || pointIsInside(center(cellId));
        if(!isBndryCell) {
          floodCells(cellId);
        }
      }
    }
  }

  void floodCells(GInt cellId) {
    std::stack<GInt> flooding;
    flooding.template emplace(cellId);
    const GBool inside = property(cellId, CellProperties::inside);

    while(!flooding.empty()) {
      const GInt currentCellId = flooding.top();
      flooding.pop();
      for(GInt id = 0; id < cartesian::maxNoNghbrs<NDIM>(); ++id) {
        const GInt nghbrId = m_nghbrIds[currentCellId].n[id];
        if(nghbrId != INVALID_CELLID && !property(nghbrId, CellProperties::marked)) {
          property(nghbrId, CellProperties::marked) = true;
          if(!property(nghbrId, CellProperties::bndry)) {
            property(nghbrId, CellProperties::inside) = inside;
            flooding.template emplace(nghbrId);
          } else {
            property(nghbrId, CellProperties::inside) = true;
          }
        }
      }
    }
  }

  [[nodiscard]] auto pointIsInside(const Point<NDIM>& x) const -> GBool { return geometry()->pointIsInside(x); }

  [[nodiscard]] auto cellHasCut(const GInt cellId) const -> GBool {
    const GDouble cellLength = lengthOnLvl(std::to_integer<GInt>(level(cellId)));
    return cellHasCut(cellId, cellLength);
  }

  [[nodiscard]] inline auto cellHasCut(const GInt cellId, const GDouble length) const -> GBool {
    return geometry()->cutWithCell(center(cellId), length);
  }

  void copyCell(const GInt from, const GInt to) {
    ASSERT(!property(from, CellProperties::markedForDeletion), "Invalid cell to be copied!");
    ASSERT(from >= 0, "Invalid from!");
    ASSERT(to >= 0, "Invalid to!");

    property(to)     = property(from);
    level(to)        = level(from);
    center(to)       = center(from);
    globalId(to)     = globalId(from);
    parent(to)       = parent(from);
    m_nghbrIds[to]   = m_nghbrIds[from];
    m_childIds[to]   = m_childIds[from];
    m_noChildren[to] = m_noChildren[from];

    for(GInt dir = 0; dir < cartesian::maxNoNghbrs<NDIM>(); ++dir) {
      if(m_nghbrIds[to].n[dir] != INVALID_CELLID) {
        m_nghbrIds[m_nghbrIds[to].n[dir]].n[cartesian::oppositeDir(dir)] = to;
      }
      m_nghbrIds[from].n[dir] = INVALID_CELLID;
    }

    if(parent(to) != INVALID_CELLID) {
      updateParent(parent(to), from, to);
    }

    for(GInt childId = 0; childId < cartesian::maxNoChildren<NDIM>(); ++childId) {
      if(m_childIds[to].c[childId] != INVALID_CELLID) {
        parent(m_childIds[to].c[childId]) = to;
      }
    }
  }

  void swapCell(GInt idA, GInt idB) {
    copyCell(idA, capacity() - 1);
    copyCell(idB, idA);
    copyCell(capacity() - 1, idB);
  }

  void updateParent(const GInt parentId, const GInt oldChildCellId, const GInt newChildCellId) {
    ASSERT(parentId >= 0, "Invalid parentId!");
    for(GInt childId = 0; childId < cartesian::maxNoChildren<NDIM>(); ++childId) {
      if(m_childIds[parentId].c[childId] == oldChildCellId) {
        m_childIds[parentId].c[childId] = newChildCellId;
        return;
      }
    }
  }

  void reorderHilbertCurve() {
    logger << SP2 << "+ reordering grid based on Hilbert curve" << std::endl;
    std::cout << SP2 << "+ reordering grid based on Hilbert curve" << std::endl;

    Point<NDIM>       centerOfGravity = Point<NDIM>(cog().data());
    std::vector<GInt> hilbertIds(size());
    std::vector<GInt> index(size());
    std::vector<GInt> pos(size());
    std::vector<GInt> rev(size());
    // generate indices
    std::iota(index.begin(), index.end(), 0);
    std::iota(pos.begin(), pos.end(), 0);
    std::iota(rev.begin(), rev.end(), 0);

    GInt hilbertLevel = partitionLvl();
    for(GInt cellId = 0; cellId < size(); ++cellId) {
      // Normalization to unit cube
      // array() since there is no scalar addition for vectors...
      Point<NDIM> x      = ((center(cellId) - centerOfGravity).array() + HALF * lengthOnLvl(0)) / lengthOnLvl(0);
      hilbertIds[cellId] = hilbert::index<NDIM>(x, hilbertLevel);
    }
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      logger << "checking duplicated Hilbert Ids" << std::endl;
      std::vector<GInt> duplicatedIds = checkDuplicateIds(hilbertIds);
      if(!duplicatedIds.empty()) {
        for(auto id : duplicatedIds) {
          GInt cellId = 0;
          std::cerr << "duplicated id " << id << std::endl;
          for(auto hilbertId : hilbertIds) {
            if(id == hilbertId) {
              std::cerr << "cellId " << cellId << std::endl;
              std::cerr << center(cellId) << std::endl;
            }
            ++cellId;
          }
        }
        TERMM(-1, "Duplicated Hilbert Ids found!");
      }
    }

    // sort index by hilbertId
    std::sort(index.begin(), index.end(), [&](int A, int B) -> bool { return hilbertIds[A] < hilbertIds[B]; });

    for(GInt id = 0; id < size(); ++id) {
      const GInt hilbertPos = index[id];
      if(id != hilbertPos) {
        // is not already in correct position due to swapping
        if(id != pos[hilbertPos]) {
          swapCell(id, pos[hilbertPos]);
        }
        GInt tmp             = rev[id];
        rev[pos[hilbertPos]] = tmp;
        pos[tmp]             = pos[hilbertPos];
      }
    }

    if(DEBUG_LEVEL > Debug_Level::debug) {
      std::fill(hilbertIds.begin(), hilbertIds.end(), 0);
      for(GInt cellId = 0; cellId < size(); ++cellId) {
        // Normalization to unit cube
        // array() since there is no scalar addition for vectors...
        Point<NDIM> x      = ((center(cellId) - centerOfGravity).array() + HALF * lengthOnLvl(0)) / lengthOnLvl(0);
        hilbertIds[cellId] = hilbert::index<NDIM>(x, hilbertLevel);
      }
      for(GInt cellId = 1; cellId < size(); ++cellId) {
        if(hilbertIds[cellId - 1] > hilbertIds[cellId]) {
          std::cerr << hilbertIds[cellId - 1] << " > " << hilbertIds[cellId] << std::endl;
          TERMM(-1, "Hilbert Ids not ordered correctly!");
        }
      }
    }
  }

  std::vector<LevelOffsetType>    m_levelOffsets{};
  std::vector<GInt>               m_noChildren{};
  std::vector<NeighborList<NDIM>> m_nghbrIds{};
  std::vector<ChildList<NDIM>>    m_childIds{};
  std::vector<GInt>               m_rfnDistance{};
};
#endif // GRIDGENERATOR_CARTESIANGRID_GENERATION_H
