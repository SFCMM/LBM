#ifndef LBM_SURFACE_H
#define LBM_SURFACE_H

#include "cartesiangrid_base.h"

// 3D: Surface 2D: Line 1D: Point
class SurfaceInterface {
 public:
  [[nodiscard]] virtual auto getCellList() const -> const std::vector<GInt>& = 0;
  virtual void               setCellList(const std::vector<GInt>& cellList)  = 0;
  virtual void               addCell(const GInt cellId, const GInt dir)      = 0;

  [[nodiscard]] virtual auto normal_p(const GInt surfCellId) const -> const GDouble*               = 0;
  [[nodiscard]] virtual auto neighbor(const GInt cellId, const GInt dir) const -> GInt             = 0;
  [[nodiscard]] virtual auto property(const GInt cellId, const CellProperties prop) const -> GBool = 0;

 private:
};

// 3D: Surface 2D: Line 1D: Point
template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class Surface : public SurfaceInterface {
 public:
  explicit Surface(CartesianGridData<NDIM> data, grid::cell::BitsetType* properties) : m_grid(data), m_properties(properties){};
  ~Surface() = default;

  Surface(const Surface<DEBUG_LEVEL, NDIM>& copy) = default;


  // todo:fix
  //  Surface(Surface&&)      = delete;
  //    auto operator=(const Surface& copy) -> Surface& {
  //      m_center = copy.m_center;
  //      m_normal = copy.m_normal;
  //      m_cellId = copy.m_cellId;
  //      m_grid = copy.m_grid;
  //    }
  //  auto operator=(Surface&&) -> Surface& = delete;


  [[nodiscard]] auto getCellList() const -> const std::vector<GInt>& override { return m_cellId; }

  void setCellList(const std::vector<GInt>& cellList) override {
    if(cellList.empty()) {
      TERMM(-1, "Invalid cellList ");
    }

    m_cellId.clear();
    //    m_nghbrIds.clear();
    std::copy(cellList.begin(), cellList.end(), std::back_inserter(m_cellId));


    //    updateNeighbors();
  }

  void addCell(const GInt cellId, const GInt dir) override {
    m_cellId.emplace_back(cellId);
    m_normal[cellId] = cartesian::dirVec<NDIM>(dir);
  }

  void removeCell(const GInt cellId) {
    m_cellId.erase(std::find(m_cellId.begin(), m_cellId.end(), cellId));
    m_normal.erase(cellId);
  }

  auto normal(const GInt surfCellId) const -> const VectorD<NDIM>& { return m_normal.at(surfCellId); }

  [[nodiscard]] auto normal_p(const GInt surfCellId) const -> const GDouble* override { return &m_normal.at(surfCellId)[0]; }


  [[nodiscard]] auto size() const -> GInt {
    return m_cellId.size();
  }

  auto center(const GInt cellId) const -> const VectorD<NDIM>& { return m_grid.center(cellId); }

  [[nodiscard]] auto cellLength(const GInt cellId) const -> GDouble { return m_grid.cellLength(cellId); }

  auto grid() const -> CartesianGridData<NDIM> { return m_grid; }

  [[nodiscard]] auto neighbor(const GInt cellId, const GInt dir) const -> GInt override {
    //    if(DEBUG_LEVEL >= Debug_Level::debug) {
    //      if(m_nghbrIds.empty()) {
    //        TERMM(-1, "Surface not initialized!");
    //      }
    //    }

    if(DEBUG_LEVEL >= Debug_Level::debug) {
      //      if(m_nghbrIds.count(cellId) == 0){
      //        cerr0<<"ERROR: Invalid cellId " << cellId << std::endl;
      //        std::exit(-1);
      //      }
      if(dir >= cartesian::maxNoNghbrsDiag<NDIM>()) {
        cerr0 << "ERROR: Invalid direction for " << cellId << " and " << dir << std::endl;
        std::exit(-1);
      }
    }
    return m_grid.neighbor(cellId, dir);
    //    return m_nghbrIds.at(cellId)[dir];
  }

  // todo: probably not needed
  void setBndryGhostCells() { m_hasBndryGhosts = true; }

  // todo: probably not needed
  [[nodiscard]] auto hasBndryGhostCells() const -> GBool { return m_hasBndryGhosts; }

  [[nodiscard]] auto property(const GInt cellId, const CellProperties prop) const -> GBool override {
    return m_grid.property(cellId, prop);
  }

  void property(const GInt cellId, const CellProperties prop, const GBool value) { m_properties[cellId][static_cast<GInt>(prop)] = value; }

  auto setProperty(const CellProperties prop, const GBool value) {
    for(const GInt surfCellId : m_cellId) {
      property(surfCellId, prop, value);
    }
  }

 private:
  std::vector<GInt> m_cellId;

  std::unordered_map<GInt, VectorD<NDIM>> m_normal;
  CartesianGridData<NDIM>                 m_grid;
  grid::cell::BitsetType*                 m_properties = nullptr;

  GBool m_hasBndryGhosts = false;
};

#endif // LBM_SURFACE_H
