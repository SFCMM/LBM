#ifndef LBM_LBM_EQUILIBRIUM_FUNC_H
#define LBM_LBM_EQUILIBRIUM_FUNC_H

#include <sfcmm_common.h>
#include "lbm_constants.h"


namespace eq {
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
