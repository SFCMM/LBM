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
    m_cellId.clear();
    std::copy(cellList.begin(), cellList.end(), std::back_inserter(m_cellId));
  }

  void addCell(const GInt cellId, const GInt dir) override {
    m_cellId.emplace_back(cellId);
    m_normal.emplace_back();
    m_normal.back().fill(0);
    m_normal.back()[dir / 2] = 2 * (dir % 2) - 1;
  }

  auto normal(const GInt surfCellId) const -> const VectorD<NDIM>& { return m_normal[surfCellId]; }

  [[nodiscard]] auto size() const -> GInt {
    ASSERT(!m_cellId.empty(), "Not inited!");
    return m_cellId.size();
  }

  auto center(const GInt cellId) const -> const VectorD<NDIM>& { return m_grid.center(cellId); }

  [[nodiscard]] auto cellLength(const GInt cellId) const -> GDouble { return m_grid.cellLength(cellId); }

  auto grid() const -> CartesianGridData<NDIM> { return m_grid; }

 private:
  std::vector<GInt>             m_cellId;
  std::vector<VectorD<NDIM>>    m_normal;
  const CartesianGridData<NDIM> m_grid = nullptr;
};

#endif // LBM_SURFACE_H
