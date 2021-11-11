#ifndef LBM_SURFACE_H
#define LBM_SURFACE_H

// 3D: Surface 2D: Line 1D: Point
class SurfaceInterface {
 public:
  [[nodiscard]] virtual auto getCellList() const -> const std::vector<GInt>& = 0;
  virtual void               setCellList(const std::vector<GInt>& cellList)  = 0;
  virtual void               addCell(const GInt cellId, const GInt dir)      = 0;

 private:
};

// 3D: Surface 2D: Line 1D: Point
template <GInt NDIM>
class Surface : public SurfaceInterface {
 public:
  Surface() = default;
  ~Surface() = default;

  //todo:fix
//  Surface(const Surface&) = delete;
//  Surface(Surface&&)      = delete;
//  auto operator=(const Surface&) -> Surface& = delete;
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
    m_normal.back()[dir / 2] = 1 - 2 * ((dir-1)%2);
  }

  auto normal(const GInt surfId) const -> const VectorD<NDIM>&{
    return m_normal[surfId];
  }

  [[nodiscard]] auto size() const -> GInt {
    ASSERT(!m_cellId.empty(), "Not inited!");
    return m_cellId.size();
  }

 private:
  std::vector<GInt>          m_cellId;
  std::vector<VectorD<NDIM>> m_normal;
};

#endif // LBM_SURFACE_H