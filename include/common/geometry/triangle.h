// SPDX-License-Identifier: BSD-3-Clause

#ifndef GRIDGENERATOR_TRIANGLE_H
#define GRIDGENERATOR_TRIANGLE_H

#include "../util/string_helper.h"

// namespace to hide local point definition
namespace triangle_ {
template <GInt NDIM>
using Point = VectorD<NDIM>;
}

template <GInt NDIM>
struct triangle {
  std::array<triangle_::Point<NDIM>, 3> m_vertices;
  triangle_::Point<NDIM>                m_normal;

  // bounding box
  triangle_::Point<NDIM> m_max;
  triangle_::Point<NDIM> m_min;
};
namespace triangle_ {
template <GInt NDIM>
[[nodiscard]] inline auto min(const triangle<NDIM>& tri, const GInt dir) -> GDouble {
  return tri.m_min[dir];
}

template <GInt NDIM>
[[nodiscard]] inline auto max(const triangle<NDIM>& tri, const GInt dir) -> GDouble {
  return tri.m_max[dir];
}

template <GInt NDIM>
[[nodiscard]] inline auto boundingBox(const triangle<NDIM>& tri, const GInt dir) -> GDouble {
  ASSERT(dir <= 2 * NDIM, "Invalid dir");

  if(dir >= NDIM) {
    return tri.m_max(dir - NDIM);
  }
  return tri.m_min(dir);
}

template <GInt NDIM>
inline void print(const triangle<NDIM>& tri) {
  std::cout << "max " << strStreamify<NDIM>(tri.m_max).str() << " min " << strStreamify<NDIM>(tri.m_min).str() << "\n";
  std::cout << "Verticies" << '\n';
  std::cout << strStreamify<NDIM>(tri.m_vertices[0]).str() << "\n";
  std::cout << strStreamify<NDIM>(tri.m_vertices[1]).str() << "\n";
  std::cout << strStreamify<NDIM>(tri.m_vertices[2]).str() << "\n";
  std::cout << "Normal"
            << "\n";
  std::cout << strStreamify<NDIM>(tri.m_normal).str() << std::endl;
}

} // namespace triangle_
#endif // GRIDGENERATOR_TRIANGLE_H
