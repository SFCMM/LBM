#ifndef LBM_MOMENTS_H
#define LBM_MOMENTS_H

#include <sfcmm_common.h>
#include "constants.h"
#include "variables.h"

template <GInt NDIM, GInt NDIST, LBEquation EQ>
inline void calcDensity(const std::vector<GInt>& cells, const std::function<GDouble&(GInt, GInt)>& fold,
                        const std::function<GDouble&(GInt, GInt)>& vars) {
  if(EQ == LBEquation::Navier_Stokes) {
    for(const GInt cellId : cells) {
      vars(cellId, LBMVariables<LBEquation::Navier_Stokes, NDIM>::rho()) = fold(cellId, 0);
      for(GInt dist = 1; dist < NDIST; ++dist) {
        vars(cellId, LBMVariables<LBEquation::Navier_Stokes, NDIM>::rho()) += fold(cellId, dist);
      }
    }
  } else {
    TERMM(-1, "Invalid LBEquation");
  }
}


#endif // LBM_MOMENTS_H
