#ifndef LBM_LBM_EQUILIBRIUM_FUNC_H
#define LBM_LBM_EQUILIBRIUM_FUNC_H

#include <sfcmm_common.h>
#include "lbm_constants.h"


namespace eq {

static inline auto defaultEq(const GDouble weight, const GDouble density, const GDouble cu, const GDouble vsq) -> GDouble {
  return weight * density * (1.0 + cu / lbm_cssq + cu * cu / (2.0 * lbm_cssq * lbm_cssq) - vsq / (2.0 * lbm_cssq));
}

template <LBMethodType LBTYPE>
static inline void defaultEq(GDouble* feq, const GDouble density, const GDouble* const velocity) {
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;

  GDouble vsq = 0;
  for(GInt dir = 0; dir < NDIM; ++dir) {
    vsq += velocity[dir] * velocity[dir];
  }

  for(GInt dist = 0; dist < NDIST; ++dist) {
    GDouble cu = 0;
    for(GInt dir = 0; dir < NDIM; ++dir) {
      cu += velocity[dir] * LBMethod<LBTYPE>::m_dirs[dist][dir];
    }
    feq[dist] = eq::defaultEq(LBMethod<LBTYPE>::m_weights[dist], density, cu, vsq);
  }
}

static inline auto poisson(const GDouble weight, const GDouble potential) -> GDouble { return weight * potential; }

template <LBMethodType LBTYPE>
static inline void poisson(GDouble* feq, const GDouble potential) {
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  for(GInt dist = 0; dist < NDIST - 1; ++dist) {
    feq[dist] = eq::poisson(LBMethod<LBTYPE>::m_weights[dist], potential);
  }
  feq[NDIST - 1] = eq::poisson((LBMethod<LBTYPE>::m_weights[NDIST - 1] - 1.0), potential);
}
} // namespace eq

#endif // LBM_LBM_EQUILIBRIUM_FUNC_H
