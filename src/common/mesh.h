#ifndef LBM_MESH_H
#define LBM_MESH_H

template <GInt NDIM>
class MeshData {
 public:
  MeshData(const GInt noMeshPoints, const std::vector<Point<NDIM>>& centers) : m_noMeshPoints(noMeshPoints), m_center(centers) {}

  [[nodiscard]] inline auto center(const GInt cellId) const -> const Point<NDIM>& { return m_center[cellId]; }

 private:
  GInt m_noMeshPoints = -1;

  const std::vector<Point<NDIM>>& m_center;
};

#endif // LBM_MESH_H
