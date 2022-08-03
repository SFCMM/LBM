#ifndef LBM_BND_WETNODE_H
#define LBM_BND_WETNODE_H

#include "common/surface.h"
#include "lbm/constants.h"
#include "lbm/moments.h"

/// Base class for wet node boundary conditions.
/// \tparam DEBUG_LEVEL
/// \tparam LBTYPE
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBnd_wallWetnode {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<LBEquationType::Navier_Stokes, NDIM>;

 public:
  LBMBnd_wallWetnode(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf) {
    // determine limited dist set etc.
    for(const GInt cellId : surf->getCellList()) {
      m_limDist.emplace_back();
      m_limConst.emplace_back();
      m_limConst.back().fill(0);
      m_normal.emplace_back(VectorD<NDIM>(surf->normal_p(cellId)));
      for(GInt dir = 0; dir < NDIST - 1; ++dir) {
        const GInt oppositeDist = LBMethod<LBTYPE>::oppositeDist(dir);

        if(inDirection<dim(LBTYPE)>(m_normal.back(), LBMethod<LBTYPE>::m_dirs[dir])
           && (surf->property(cellId, CellProperties::periodic) || surf->neighbor(cellId, oppositeDist) != INVALID_CELLID)) {
          m_limDist.back().emplace(dir);

          m_limConst.back()[dir] = 2;

        } else if(orthogonal<dim(LBTYPE)>(m_normal.back(), LBMethod<LBTYPE>::m_dirs[dir])
                  && (surf->property(cellId, CellProperties::periodic) || surf->neighbor(cellId, oppositeDist) != INVALID_CELLID)) {
          m_limDist.back().emplace(dir);

          // cell is a corner
          if(surf->neighbor(cellId, dir) == INVALID_CELLID && !surf->property(cellId, CellProperties::periodic)) {
            // this is only correct for wall velocity = 0
            // and if it is not a periodic cell will otherwise lead to reflections
            m_limConst.back()[dir] = 2;
          } else {
            m_limConst.back()[dir] = 1;
          }
        }
      }
      const GInt sumC = std::accumulate(m_limConst.back().begin(), m_limConst.back().end(), 0);
      if(sumC != NDIST - 1) {
        // we clear to mark a corner
        m_limDist.back().clear();
      } else {
        m_limDist.back().emplace(NDIST - 1);
        m_limConst.back()[NDIST - 1] = 1;
      }
    }
  }

  [[nodiscard]] auto limitedDist() const -> const std::vector<std::set<GInt>>& { return m_limDist; }

  auto limitedConst() const -> const std::vector<std::array<GDouble, NDIST>>& { return m_limConst; }

  auto normal() const -> const std::vector<VectorD<NDIM>>& { return m_normal; }

 private:
  std::vector<std::set<GInt>>             m_limDist;
  std::vector<std::array<GDouble, NDIST>> m_limConst;
  std::vector<VectorD<NDIM>>              m_normal;
};

#endif // LBM_BND_WETNODE_H
