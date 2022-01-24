#ifndef LBM_SURFACE_H
#define LBM_SURFACE_H

#include "cartesiangrid_base.h"

// 3D: Surface 2D: Line 1D: Point
class SurfaceInterface {
 public:
  [[nodiscard]] virtual auto getCellList() const -> const std::vector<GInt>& = 0;
  virtual void               setCellList(const std::vector<GInt>& cellList)  = 0;
  virtual void               addCell(const GInt cellId, const GInt dir)      = 0;

 private:
};

// 3D: Surface 2D: Line 1D: Point
template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class Surface : public SurfaceInterface {
 public:
  Surface(CartesianGridData<NDIM> data) : m_grid(data){};
  ~Surface() = default;

  // todo:fix
  //  Surface(const Surface& copy) : m_center(copy.m_center), m_grid(copy.m_grid) {}
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
    m_nghbrIds.clear();
    std::copy(cellList.begin(), cellList.end(), std::back_inserter(m_cellId));


    updateNeighbors();
  }

  void addCell(const GInt cellId, const GInt dir) override {
    m_cellId.emplace_back(cellId);
    m_normal[cellId] = cartesian::dirVec<NDIM>(dir);
  }

  void removeCell(const GInt cellId) {
    m_cellId.erase(std::find(m_cellId.begin(), m_cellId.end(), cellId));
    m_normal.erase(cellId);
  }

  void updateNeighbors() {
    for(const auto cellId : m_cellId) {
      std::array<GInt, cartesian::maxNoNghbrsDiag<NDIM>()> tmpNghbr;
      for(GInt nghbrId = 0; nghbrId < cartesian::maxNoNghbrsDiag<NDIM>(); ++nghbrId) {
        tmpNghbr[nghbrId] = m_grid.neighbor(cellId, nghbrId);
      }
      m_nghbrIds.insert({cellId, tmpNghbr});
    }
  }

  auto normal(const GInt surfCellId) const -> const VectorD<NDIM>& { return m_normal.at(surfCellId); }

  [[nodiscard]] auto size() const -> GInt {
    ASSERT(!m_cellId.empty(), "Not inited!");
    return m_cellId.size();
  }

  auto center(const GInt cellId) const -> const VectorD<NDIM>& { return m_grid.center(cellId); }

  [[nodiscard]] auto cellLength(const GInt cellId) const -> GDouble { return m_grid.cellLength(cellId); }

  auto grid() const -> CartesianGridData<NDIM> { return m_grid; }

  [[nodiscard]] auto neighbor(const GInt cellId, const GInt dir) const -> GInt {
    if(DEBUG_LEVEL >= Debug_Level::debug) {
      if(m_nghbrIds.empty()) {
        TERMM(-1, "Surface not initialized!");
      }
    }

    return m_nghbrIds.at(cellId)[dir];
  }

  void setBndryGhostCells() { m_hasBndryGhosts = true; }

  [[nodiscard]] auto hasBndryGhostCells() const -> GBool { return m_hasBndryGhosts; }

 private:
  std::vector<GInt>                                                              m_cellId;
  std::unordered_map<GInt, std::array<GInt, cartesian::maxNoNghbrsDiag<NDIM>()>> m_nghbrIds;

  std::unordered_map<GInt, VectorD<NDIM>> m_normal;
  const CartesianGridData<NDIM>           m_grid = nullptr;

  GBool m_hasBndryGhosts = false;
};

#endif // LBM_SURFACE_H
