#ifndef LBM_LINE_H
#define LBM_LINE_H
#include <Eigen/Geometry>
#include <sfcmm_common.h>


template <GInt NDIM>
class Line {
 private:
  using EigenLine = Eigen::Hyperplane<GDouble, NDIM>;

 public:
  Line(const Point<NDIM>& pointA, const Point<NDIM>& pointB) : m_a(pointA), m_b(pointB), m_line(EigenLine::Through(pointA, pointB)) {}

  auto distance(const Point<NDIM>& point) const -> GDouble { return m_line.absDistance(point); }

 private:
  Point<NDIM> m_a;
  Point<NDIM> m_b;
  EigenLine   m_line;
};

#endif // LBM_LINE_H
