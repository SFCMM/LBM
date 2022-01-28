#ifndef LBM_MOMENTS_H
#define LBM_MOMENTS_H

#include <array>
#include <set>
#include <sfcmm_common.h>
#include "constants.h"
#include "variables.h"

/// Calculate the density for a given set of cells.
/// \tparam NDIM Dimensionality
/// \tparam NDIST Number of distributions
/// \tparam EQ equation
/// \param cells List of cells for which to calculate the density
/// \param fold Distribution to use
/// \param vars Variable storage
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

/// Calculate the density for a cell using a reduced set of distribution e.g. at boundaries
/// \tparam NDIM Dimensionality
/// \tparam NDIST Number of distributions
/// \tparam EQ Equation
/// \tparam SLIP NOSLIP velocity = 0
/// \param cellId CellId for which to calculate the density
/// \param limitedDists The available dists
/// \param constants Constants to solve the limited set.
/// \param normal Normal direction of the missing distributions
/// \param fold Distributions to use
/// \param vars Variable storage
template <GInt NDIM, GInt NDIST, LBEquation EQ, GBool NOSLIP>
inline void calcDensity_limited(const GInt cellId, const std::set<GInt>& limitedDists, const std::array<GDouble, NDIST>& constants,
                                const GDouble* normal, const std::function<GDouble&(GInt, GInt)>& fold,
                                const std::function<GDouble&(GInt, GInt)>& vars) {
  if(EQ == LBEquation::Navier_Stokes) {
    vars(cellId, LBMVariables<LBEquation::Navier_Stokes, NDIM>::rho()) = 0;
    for(const GInt dist : limitedDists) {
      vars(cellId, LBMVariables<LBEquation::Navier_Stokes, NDIM>::rho()) += constants[dist] * fold(cellId, dist);
    }

    // for no slip walls the velocity in the cell is 0!
    if(!NOSLIP) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        if(normal[dir] > GDoubleEps) {
          vars(cellId, LBMVariables<LBEquation::Navier_Stokes, NDIM>::rho()) *=
              1.0 / (1.0 + vars(cellId, LBMVariables<LBEquation::Navier_Stokes, NDIM>::velocity(dir)));
        } else if(normal[dir] < 0) {
          vars(cellId, LBMVariables<LBEquation::Navier_Stokes, NDIM>::rho()) *=
              1.0 / (1.0 - vars(cellId, LBMVariables<LBEquation::Navier_Stokes, NDIM>::velocity(dir)));
        }
      }
    }
  } else {
    TERMM(-1, "Invalid LBEquation");
  }
}


#endif // LBM_MOMENTS_H
