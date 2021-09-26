#ifndef GRIDGENERATOR_CARTESIANGRID_GENERATION_H
#define GRIDGENERATOR_CARTESIANGRID_GENERATION_H
#include "base_cartesiangrid.h"

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

  using PropertyBitsetType = grid::cell::BitsetType;
  using ChildListType      = std::array<GInt, cartesian::maxNoChildren<NDIM>()>;

  explicit CartesianGridGen(const GInt maxNoCells) : m_capacity(maxNoCells) { setCapacity(maxNoCells); }
  ~CartesianGridGen() override              = default;
  CartesianGridGen(const CartesianGridGen&) = delete;
  CartesianGridGen(CartesianGridGen&&)      = delete;
  auto operator=(const CartesianGridGen&) -> CartesianGridGen& = delete;
  auto operator=(CartesianGridGen&&) -> CartesianGridGen& = delete;

  void setCapacity(const GInt capacity) override {
    if(!m_center.empty()) {
      TERMM(-1, "Invalid operation tree already allocated.");
    }
    if(capacity < m_capacity) {
      logger << "WARNING: memory allocation requested. But capacity is sufficient CartesianGridGen::setCapacity() requested " << capacity
             << " but allocated " << m_capacity << std::endl;
      return;
    }
    m_center.resize(capacity);
    m_parentId.resize(capacity);
    m_globalId.resize(capacity);
    m_noChildren.resize(capacity);
    m_nghbrIds.resize(capacity);
    m_childIds.resize(capacity);
    m_rfnDistance.resize(capacity);
    m_properties.resize(capacity);
    m_level.resize(capacity);
    m_capacity = capacity;
  }

  void reset() override {
    m_center.clear();
    m_parentId.clear();
    m_globalId.clear();
    m_noChildren.clear();
    m_nghbrIds.clear();
    m_childIds.clear();
    m_rfnDistance.clear();
    m_properties.clear();
    m_level.clear();
  }

  void setMaxLvl(const GInt _maxLvl) override {
    m_levelOffsets.resize(_maxLvl + 1);
    BaseCartesianGrid<DEBUG_LEVEL, NDIM>::setMaxLvl(_maxLvl);
  }

  void createPartitioningGrid(const GInt partitioningLvl) override {
    RECORD_TIMER_START(TimeKeeper[Timers::GridPart]);
    if(m_capacity < 1) {
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
      m_levelOffsets[1] = {m_capacity - cartesian::maxNoChildren<NDIM>(), m_capacity};
    } else {
      // initial cell placed at the end
      m_levelOffsets[0] = {m_capacity - 1, m_capacity};
      m_levelOffsets[1] = {0, cartesian::maxNoChildren<NDIM>()};
    }

    const GInt begin                       = m_levelOffsets[0].begin;
    m_center[begin]                        = Point<NDIM>(cog().data());
    m_globalId[begin]                      = begin;
    property(begin, CellProperties::bndry) = true;
    m_size                                 = 1;

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
        const GInt newLevelBegin = m_capacity - (prevLevelNoCells)*cartesian::maxNoChildren<NDIM>();
        m_levelOffsets[l + 1]    = {newLevelBegin, m_capacity};

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

    std::fill(m_parentId.begin(), m_parentId.end(), INVALID_CELLID);
    reorderHilbertCurve();

    RECORD_TIMER_STOP(TimeKeeper[Timers::GridPart]);
  }

  void uniformRefineGrid(const GInt uniformLevel) override {
    RECORD_TIMER_START(TimeKeeper[Timers::GridUniform]);
    logger << SP1 << "Uniformly refine grid to level " << uniformLevel << std::endl;
    std::cout << SP1 << "Uniformly refine grid to level " << uniformLevel << std::endl;

    if(partitionLvl() == uniformLevel) {
      return;
    }

    for(GInt l = partitionLvl(); l < uniformLevel; l++) {
      m_levelOffsets[l + 1] = {m_size, m_size + levelSize(m_levelOffsets[l]) * cartesian::maxNoChildren<NDIM>()};
      if(m_levelOffsets[l + 1].end > m_capacity) {
        outOfMemory(l + 1);
      }

      refineGrid<true>(l);
    }
    RECORD_TIMER_STOP(TimeKeeper[Timers::GridUniform]);
  }

  void refineMarkedCells(const GInt noCellsToRefine) override {
    if(noCellsToRefine == 0) {
      logger << "WARNING: refineMarkedCells called but nothing to do" << std::endl;
      return;
    }

    logger << SP1 << "Refining marked cells to level " << currentHighestLvl() + 1 << std::endl;
    std::cout << SP1 << "Refining marked cells to level " << currentHighestLvl() + 1 << std::endl;
    // update the offsets
    m_levelOffsets[currentHighestLvl() + 1] = {m_size, m_size + noCellsToRefine * cartesian::maxNoChildren<NDIM>()};
    if(m_levelOffsets[currentHighestLvl() + 1].end > m_capacity) {
      outOfMemory(currentHighestLvl() + 1);
    }
    logger << SP2 << "* cells to refine: " << noCellsToRefine << std::endl;
    std::cout << SP2 << "* cells to refine: " << noCellsToRefine << std::endl;


    refineGrid<true, false>(currentHighestLvl());
  }

  auto markBndryCells() -> GInt override {
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

  void save(const GString& fileName, const json& gridOutConfig) override {
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
    std::function<GBool(GInt)> isLowestLevel = [&](GInt cellId) { return std::to_integer<GInt>(m_level[cellId]) == partitionLvl(); };

    // only output the highest level
    std::function<GBool(GInt)> isHighestLevel = [&](GInt cellId) { return std::to_integer<GInt>(m_level[cellId]) == currentHighestLvl(); };

    // only output leaf cells (i.e. cells without children)
    std::function<GBool(GInt)> isLeaf = [&](GInt cellId) { return m_noChildren[cellId] == 0; };

    // only output the lowest level
    std::function<GBool(GInt)> isTargetLevel = [&](GInt cellId) { return std::to_integer<GInt>(m_level[cellId]) == outputLvl; };


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
    } else {
      TERMM(-1, "Unknown output filter " + filter);
    }

    std::vector<GString>              index;
    std::vector<std::vector<GString>> values;
    cerr0 << "Selected output values:";
    for(const auto& outputvalue : outvalues) {
      if(outputvalue == "level") {
        index.emplace_back("Level");
        values.emplace_back(toStringVector(m_level, m_size));
        cerr0 << " level ";
      } else if(outputvalue == "noChildren") {
        index.emplace_back("NoChildren");
        values.emplace_back(toStringVector(m_noChildren, m_size));
        cerr0 << " noChildren ";
      } else {
        logger << "WARNING: The output value " + outputvalue + " is not a valid output!" << std::endl;
      }
    }
    cerr0 << std::endl;

    if(format == "ASCII") {
      ASCII::writePointsCSV<NDIM>(fileName, m_size, m_center, index, values, outputFilter);
    } else if(format == "VTK") {
      VTK::ASCII::writePoints<NDIM>(fileName, m_size, m_center, index, values, outputFilter);
    } else if(format == "VTKB") {
      // todo: rename format
      VTK::BINARY::writePoints<NDIM>(fileName, m_size, m_center, index, values, outputFilter);
    } else {
      TERMM(-1, "Unknown output format " + format);
    }
  }

  static constexpr auto memorySizePerCell() -> GInt {
    return sizeof(GInt) * (1 + 1 + 1 + 1 + 2) // m_parentId, m_globalId, m_noChildren, m_rfnDistance, m_levelOffsets
           + sizeof(Point<NDIM>)              // m_center
           + sizeof(NeighborList<NDIM>)       // m_nghbrIds
           + sizeof(ChildList<NDIM>)          // m_childIds
           + sizeof(PropertyBitsetType)       // m_properties
           + 1;                               // m_level
  }


 private:
  inline auto                  property(const GInt id, CellProperties p) -> auto { return m_properties[id][static_cast<GInt>(p)]; }
  [[nodiscard]] inline auto    property(const GInt id, CellProperties p) const -> GBool { return m_properties[id][static_cast<GInt>(p)]; }
  std::vector<LevelOffsetType> m_levelOffsets{};
  std::vector<Point<NDIM>>     m_center{};
  std::vector<GInt>            m_parentId{INVALID_CELLID};
  std::vector<GInt>            m_globalId{INVALID_CELLID};
  std::vector<GInt>            m_noChildren{};
  std::vector<std::byte>       m_level{};
  std::vector<NeighborList<NDIM>> m_nghbrIds{};
  std::vector<ChildList<NDIM>>    m_childIds{};
  std::vector<GInt>               m_rfnDistance{};
  std::vector<PropertyBitsetType> m_properties{};

  GInt m_capacity{0};
  GInt m_size{0};

  void outOfMemory(GInt level) {
    cerr0 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    logger << "ERROR: Not enough memory to generate grid! Increase maxNoCells: " << m_capacity << std::endl;
    cerr0 << "ERROR: Not enough memory to generate grid! Increase maxNoCells: " << m_capacity << std::endl;
    logger << "level " << level - 1 << " [" << m_levelOffsets[level - 1].begin << ", " << m_levelOffsets[level - 1].end << "]" << std::endl;
    cerr0 << "level " << level - 1 << " [" << m_levelOffsets[level - 1].begin << ", " << m_levelOffsets[level - 1].end << "]" << std::endl;
    logger << "level " << level << " [" << m_levelOffsets[level].begin << ", " << m_levelOffsets[level].end << "]" << std::endl;
    cerr0 << "level " << level << " [" << m_levelOffsets[level].begin << ", " << m_levelOffsets[level].end << "]" << std::endl;
    cerr0 << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

    TERMM(-1, "Out of memory!");
  }

  template <GBool DISTRIBUTED = false, GBool UNIFORM = true>
  void refineGrid(const GInt level) {
    if(DISTRIBUTED && !MPI::isSerial()) {
      // updateHaloOffsets();
    }

    if(UNIFORM) {
      refineGrid(m_levelOffsets, level);
      if(DISTRIBUTED && !MPI::isSerial()) {
        // todo:implement
        // refineGrid(m_haloOffsets, l);
      }
    } else {
      refineGridMarkedOnly(m_levelOffsets, level);
      if(DISTRIBUTED && !MPI::isSerial()) {
        // todo:implement
        // refineGridMarkedOnly(m_haloOffsets, l);
      }
    }
    m_size = m_levelOffsets[level + 1].end;

    findChildLevelNghbrs(m_levelOffsets, level);
    if(DISTRIBUTED && !MPI::isSerial()) {
      // todo:implement
      // findChildLevelNeighbors(m_haloOffsets, l);
    }

    if(DISTRIBUTED && !MPI::isSerial()) {
      // todo:implement
      //        deleteOutsideCellsParallel(l + 1);
    } else {
      deleteOutsideCells(level + 1);
    }

    m_size = m_levelOffsets[level + 1].end;
    increaseCurrentHighestLvl();
    logger.updateAttributes();
  }

  void refineGrid(const std::vector<LevelOffsetType>& levelOffset, const GInt level) {
    logger << SP2 << "+ refining grid on level: " << level << std::endl;
    std::cout << SP2 << "+ refining grid on level: " << level << std::endl;

    ASSERT(level + 1 <= maxLvl(), "Invalid refinement level! " + std::to_string(level + 1) + ">" + std::to_string(maxLvl()));

    // refine all cells on the given level
#ifdef _OPENMP
#pragma omp parallel default(none) shared(levelOffset, level)
    {
#endif

      const GInt begin = levelOffset[level].begin;
      const GInt end   = levelOffset[level].end;
#ifdef _OPENMP
#pragma omp for
#endif
      for(GInt cellId = begin; cellId < end; ++cellId) {
        const GInt cellCount = cellId - begin;
        refineCell(cellId, levelOffset[level + 1].begin + cellCount * cartesian::maxNoChildren<NDIM>());
      }
#ifdef _OPENMP
    }
#endif
  }

  void refineGridMarkedOnly(const std::vector<LevelOffsetType>& levelOffset, const GInt level) {
    GInt refinedCells = 0;
    for(GInt cellId = levelOffset[level].begin; cellId < levelOffset[level].end; ++cellId) {
      if(property(cellId, CellProperties::toRefine)) {
        refineCell(cellId, levelOffset[level + 1].begin + refinedCells * cartesian::maxNoChildren<NDIM>());
        ++refinedCells;
      }
    }
  }

  void refineCell(const GInt cellId, const GInt offset) {
    if(DEBUG_LEVEL > Debug_Level::debug) {
      logger << "refine cell " << cellId << " with offset " << offset << std::endl;
    }
    const GInt    refinedLvl       = std::to_integer<GInt>(m_level[cellId]) + 1;
    const GDouble refinedLvlLength = lengthOnLvl(refinedLvl);

    for(GInt childId = 0; childId < cartesian::maxNoChildren<NDIM>(); ++childId) {
      const GInt childCellId = offset + childId;
      m_center[childCellId] =
          m_center[cellId]
          + HALF * refinedLvlLength
                * Point<NDIM>(cartesian::childDir[childId].data()); // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
      m_level[childCellId]    = static_cast<std::byte>(refinedLvl);
      m_parentId[childCellId] = cellId;
      m_globalId[childCellId] = childCellId;

      // reset since we overwrite previous levels
      m_noChildren[childCellId] = 0;
      m_properties[childCellId].reset();
      m_childIds[childCellId] = {INVALID_LIST<cartesian::maxNoChildren<NDIM>()>()};
      m_nghbrIds[childCellId] = {INVALID_LIST<cartesian::maxNoNghbrs<NDIM>()>()};

      // if parent is a boundary cell check for children as well
      if(property(cellId, CellProperties::bndry)) {
        property(childCellId, CellProperties::bndry) = cellHasCut(childCellId);
      }

      // update parent
      m_childIds[cellId].c[childId] = childCellId;
      m_noChildren[cellId]++;
    }
  }

  void findChildLevelNghbrs(const std::vector<LevelOffsetType>& levelOffset, const GInt level) {
    // check all children at the given level
    for(GInt parentId = levelOffset[level].begin; parentId < levelOffset[level].end; ++parentId) {
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

            if(DEBUG_LEVEL > Debug_Level::min_debug && neighbors[dir] != INVALID_CELLID
               && (m_center[neighbors[dir]] - m_center[children[childId]]).norm()
                      > 1.9 * lengthOnLvl(std::to_integer<GInt>(m_level[children[childId]]))) {
              cerr0 << "neighbors[dir] " << neighbors[dir] << " cellId " << children[childId] << std::endl;
              cerr0 << "neighbors " << strStreamify<NDIM>(m_center[neighbors[dir]]).str() << std::endl;
              cerr0 << "neighbors " << strStreamify<NDIM>(m_center[children[childId]]).str() << std::endl;
              cerr0 << "ndiff " << (m_center[neighbors[dir]] - m_center[children[childId]]).norm() << " vs "
                    << lengthOnLvl(std::to_integer<GInt>(m_level[children[childId]])) << std::endl;
              cerr0 << "parentId " << parentId << " np " << m_nghbrIds[parentId].n[dir] << std::endl;
              cerr0 << "parent " << strStreamify<NDIM>(m_center[parentId]).str() << std::endl;
              cerr0 << "parent neighbors " << strStreamify<NDIM>(m_center[m_nghbrIds[parentId].n[dir]]).str() << std::endl;
              cerr0 << "pdiff " << (m_center[parentId] - m_center[m_nghbrIds[parentId].n[dir]]).norm() << " vs "
                    << lengthOnLvl(std::to_integer<GInt>(m_level[parentId])) << std::endl;

              TERMM(-1, "Invalid neighbor");
            }
          }
        }
      }
    }
  }

  void deleteOutsideCells(const GInt level) {
    markOutsideCells(m_levelOffsets, level);

    // delete cells that have been marked as being outside
    for(GInt cellId = m_levelOffsets[level].end - 1; cellId >= m_levelOffsets[level].begin; --cellId) {
      ASSERT(!property(cellId, CellProperties::bndry)
                 || property(cellId, CellProperties::inside) == property(cellId, CellProperties::bndry),
             "Properties not set correctly! bndry implies IsInside!");
      // remove cell since it is not inside
      if(!property(cellId, CellProperties::inside)) {
        const GInt parentId = m_parentId[cellId];
        ASSERT(parentId != INVALID_CELLID, "Invalid parentId! (cellId: " + std::to_string(cellId) + ")");

        // remove from parent
        updateParent(m_parentId[cellId], cellId, INVALID_CELLID);
        --m_noChildren[parentId];

        // remove from neighbors
        for(GInt dir = 0; dir < cartesian::maxNoNghbrs<NDIM>(); ++dir) {
          const GInt nghbrCellId = m_nghbrIds[cellId].n[dir];
          if(nghbrCellId != INVALID_CELLID) {
            m_nghbrIds[nghbrCellId].n[cartesian::oppositeDir(dir)] = INVALID_CELLID;
          }
        }
        if(cellId != m_levelOffsets[level].end - 1) {
          // copy an inside cell to the current position to fill the hole
          copyCell(m_levelOffsets[level].end - 1, cellId);
        }
        m_levelOffsets[level].end--;
      }
    }
    m_size = levelSize(m_levelOffsets[level]);
    logger << SP3 << "* grid has " << m_size << " cells" << std::endl;
    std::cout << SP3 << "* grid has " << m_size << " cells" << std::endl;
  }

  void markOutsideCells(const std::vector<LevelOffsetType>& levelOffset, const GInt level) {
    for(GInt cellId = levelOffset[level].begin; cellId < levelOffset[level].end; ++cellId) {
      property(cellId, CellProperties::marked) = false;
    }

    for(GInt cellId = levelOffset[level].begin; cellId < levelOffset[level].end; ++cellId) {
      if(property(cellId, CellProperties::marked)) {
        continue;
      }
      property(cellId, CellProperties::marked) = true;
      const GBool isBndryCell                  = property(cellId, CellProperties::bndry);
      property(cellId, CellProperties::inside) = isBndryCell || pointIsInside(m_center[cellId]);
      if(!isBndryCell) {
        floodCells(cellId);
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

  [[nodiscard]] auto cellHasCut(GInt cellId) const -> GBool {
    const GDouble cellLength = lengthOnLvl(std::to_integer<GInt>(m_level[cellId]));
    return geometry()->cutWithCell(m_center[cellId], cellLength);
  }

  void copyCell(const GInt from, const GInt to) {
    ASSERT(!property(from, CellProperties::markedForDeletion), "Invalid cell to be copied!");
    ASSERT(from >= 0, "Invalid from!");
    ASSERT(to >= 0, "Invalid to!");

    m_properties[to] = m_properties[from];
    m_level[to]      = m_level[from];
    m_center[to]     = m_center[from];
    m_globalId[to]   = m_globalId[from];
    m_parentId[to]   = m_parentId[from];
    m_nghbrIds[to]   = m_nghbrIds[from];
    m_childIds[to]   = m_childIds[from];
    m_noChildren[to] = m_noChildren[from];

    for(GInt dir = 0; dir < cartesian::maxNoNghbrs<NDIM>(); ++dir) {
      if(m_nghbrIds[to].n[dir] != INVALID_CELLID) {
        m_nghbrIds[m_nghbrIds[to].n[dir]].n[cartesian::oppositeDir(dir)] = to;
      }
      m_nghbrIds[from].n[dir] = INVALID_CELLID;
    }

    if(m_parentId[to] != INVALID_CELLID) {
      updateParent(m_parentId[to], from, to);
    }

    for(GInt childId = 0; childId < cartesian::maxNoChildren<NDIM>(); ++childId) {
      if(m_childIds[to].c[childId] != INVALID_CELLID) {
        m_parentId[m_childIds[to].c[childId]] = to;
      }
    }
  }

  void swapCell(GInt idA, GInt idB) {
    copyCell(idA, m_capacity - 1);
    copyCell(idB, idA);
    copyCell(m_capacity - 1, idB);
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
    std::vector<GInt> hilbertIds(m_size);
    std::vector<GInt> index(m_size);
    std::vector<GInt> pos(m_size);
    std::vector<GInt> rev(m_size);
    // generate indices
    std::iota(index.begin(), index.end(), 0);
    std::iota(pos.begin(), pos.end(), 0);
    std::iota(rev.begin(), rev.end(), 0);

    GInt hilbertLevel = partitionLvl();
    for(GInt cellId = 0; cellId < m_size; ++cellId) {
      // Normalization to unit cube
      // array() since there is no scalar addition for vectors...
      Point<NDIM> x      = ((m_center[cellId] - centerOfGravity).array() + HALF * lengthOnLvl(0)) / lengthOnLvl(0);
      hilbertIds[cellId] = hilbert::index<NDIM>(x, hilbertLevel);
    }
    if(DEBUG_LEVEL > Debug_Level::min_debug) {
      logger << "checking duplicated Hilbert Ids" << std::endl;
      std::vector<GInt> duplicatedIds = checkDuplicateIds(hilbertIds);
      if(!duplicatedIds.empty()) {
        for(auto id : duplicatedIds) {
          GInt cellId = 0;
          std::cerr << "duplicated id " << id << std::endl;
          for(auto hilbertId : hilbertIds) {
            if(id == hilbertId) {
              std::cerr << "cellId " << cellId << std::endl;
              std::cerr << m_center[cellId] << std::endl;
            }
            ++cellId;
          }
        }
        TERMM(-1, "Duplicated Hilbert Ids found!");
      }
    }

    // sort index by hilberId
    std::sort(index.begin(), index.end(), [&](int A, int B) -> bool { return hilbertIds[A] < hilbertIds[B]; });

    for(GInt id = 0; id < m_size; ++id) {
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
      for(GInt cellId = 0; cellId < m_size; ++cellId) {
        // Normalization to unit cube
        // array() since there is no scalar addition for vectors...
        Point<NDIM> x      = ((m_center[cellId] - centerOfGravity).array() + HALF * lengthOnLvl(0)) / lengthOnLvl(0);
        hilbertIds[cellId] = hilbert::index<NDIM>(x, hilbertLevel);
      }
      for(GInt cellId = 1; cellId < m_size; ++cellId) {
        if(hilbertIds[cellId - 1] > hilbertIds[cellId]) {
          std::cerr << hilbertIds[cellId - 1] << " > " << hilbertIds[cellId] << std::endl;
          TERMM(-1, "Hilbert Ids not ordered correctly!");
        }
      }
    }
  }
};
#endif // GRIDGENERATOR_CARTESIANGRID_GENERATION_H
