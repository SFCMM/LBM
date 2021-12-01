#ifndef LBM_LINE_H
#define LBM_LINE_H
#include <Eigen/Geometry>


template<GInt NDIM>
class Line{
 private:
  using EigenLine = Eigen::Hyperplane<GDouble, 2>;
 public:
  Line(const Point<NDIM>& A, const Point<NDIM>& B): m_a(A), m_b(B){
    m_line = EigenLine::Through(A, B);
  }

  auto distance(const Point<NDIM>& point) const -> GDouble{
    return m_line.absDistance(point);
  }

 private:
  Point<NDIM> m_a;
  Point<NDIM> m_b;
  EigenLine m_line;
};

#endif // LBM_LINE_H
