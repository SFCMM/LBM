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
template <GInt NDIM, GInt NDIST, LBEquationType EQ>
inline void calcVelocity(const std::vector<GInt>& cells, const std::function<GDouble&(GInt, GInt)>& fold,
                         const std::function<GDouble&(GInt, GInt)>& vars) {
  using METH = LBMethod<getLBMethodType(NDIM, NDIST)>;

  if(EQ == LBEquationType::Navier_Stokes) {
    for(const GInt cellId : cells) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(cellId, LBMVariables<EQ, NDIM>::velocity(dir)) = 0;
        for(GInt dist = 0; dist < NDIST - 1; ++dist) {
          vars(cellId, LBMVariables<EQ, NDIM>::velocity(dir)) += METH::m_dirs[dist][dir] * fold(cellId, dist);
        }
        vars(cellId, LBMVariables<EQ, NDIM>::velocity(dir)) /= vars(cellId, LBMVariables<EQ, NDIM>::rho());
      }
    }
  } else {
    TERMM(-1, "Invalid LBEquation");
  }
}

/// Calculate the density for a given set of cells.
/// \tparam NDIM Dimensionality
/// \tparam NDIST Number of distributions
/// \tparam EQ equation
/// \param cells List of cells for which to calculate the density
/// \param fold Distribution to use
/// \param vars Variable storage
template <GInt NDIM, GInt NDIST, LBEquationType EQ>
inline auto calcVelocity(const GInt cellId, const std::function<GDouble&(GInt, GInt)>& fold,
                         const std::function<GDouble&(GInt, GInt)>& vars) -> VectorD<NDIM> {
  using METH = LBMethod<getLBMethodType(NDIM, NDIST)>;

  VectorD<NDIM> tempV;
  tempV.fill(0);
  if(EQ == LBEquationType::Navier_Stokes) {
    for(GInt dir = 0; dir < NDIM; ++dir) {
      for(GInt dist = 0; dist < NDIST - 1; ++dist) {
        tempV[dir] += METH::m_dirs[dist][dir] * fold(cellId, dist);
      }
      tempV[dir] /= vars(cellId, LBMVariables<EQ, NDIM>::rho());
    }
    return tempV;
  }
  TERMM(-1, "Invalid LBEquation");
}

/// Calculate the density for a given set of cells.
/// \tparam NDIM Dimensionality
/// \tparam NDIST Number of distributions
/// \tparam EQ equation
/// \param cellId List of cells for which to calculate the density
/// \param fold Distribution to use
/// \param vars Variable storage
template <GInt NDIM, GInt NDIST, LBEquationType EQ>
inline void calcDensity(const GInt cellId, const std::function<GDouble&(GInt, GInt)>& fold,
                        const std::function<GDouble&(GInt, GInt)>& vars) {
  if(EQ == LBEquationType::Navier_Stokes) {
    vars(cellId, LBMVariables<EQ, NDIM>::rho()) = fold(cellId, 0);
    for(GInt dist = 1; dist < NDIST; ++dist) {
      vars(cellId, LBMVariables<EQ, NDIM>::rho()) += fold(cellId, dist);
    }

  } else if(EQ == LBEquationType::Poisson) {
    using METH = LBMethod<getLBMethodType(NDIM, NDIST)>;

    GDouble acc_f = 0;
    for(GInt dist = 0; dist < NDIST - 1; ++dist) {
      acc_f += fold(cellId, dist);
    }
    vars(cellId, LBMVariables<EQ, NDIM>::electricPotential()) = 1.0 / (1.0 - METH::m_weights[NDIST - 1]) * acc_f;

  } else {
    TERMM(-1, "Invalid LBEquation");
  }
}

/// Calculate the density for a given set of cells.
/// \tparam NDIM Dimensionality
/// \tparam NDIST Number of distributions
/// \tparam EQ equation
/// \param cells List of cells for which to calculate the density
/// \param fold Distribution to use
/// \param vars Variable storage
template <GInt NDIM, GInt NDIST, LBEquationType EQ>
inline void calcDensity(const std::vector<GInt>& cells, const std::function<GDouble&(GInt, GInt)>& fold,
                        const std::function<GDouble&(GInt, GInt)>& vars) {
  for(const GInt cellId : cells) {
    calcDensity<NDIM, NDIST, EQ>(cellId, fold, vars);
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
template <GInt NDIM, GInt NDIST, LBEquationType EQ, GBool NOSLIP>
inline void calcDensity_limited(const GInt cellId, const std::set<GInt>& limitedDists, const std::array<GDouble, NDIST>& constants,
                                const GDouble* normal, const std::function<GDouble&(GInt, GInt)>& fold,
                                const std::function<GDouble&(GInt, GInt)>& vars) {
  if(limitedDists.empty()) {
    // skipping empty distributions which for example occur in corners of non-periodic bnd cells
    return;
  }

  if(EQ == LBEquationType::Navier_Stokes) {
    //    const GInt sumC = std::accumulate(constants.begin(), constants.end(), 0);

    //    const GDouble tmpRho = vars(cellId, LBMVariables<LBEquationType::Navier_Stokes, NDIM>::rho());
    vars(cellId, LBMVariables<LBEquationType::Navier_Stokes, NDIM>::rho()) = 0;
    for(const GInt dist : limitedDists) {
      vars(cellId, LBMVariables<LBEquationType::Navier_Stokes, NDIM>::rho()) += constants[dist] * fold(cellId, dist);
    }

    // for no slip walls the velocity in the cell is 0!
    if(!NOSLIP) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        if(normal[dir] > GDoubleEps) {
          vars(cellId, LBMVariables<LBEquationType::Navier_Stokes, NDIM>::rho()) *=
              1.0 / (1.0 + vars(cellId, LBMVariables<LBEquationType::Navier_Stokes, NDIM>::velocity(dir)));
        } else if(normal[dir] < 0) {
          vars(cellId, LBMVariables<LBEquationType::Navier_Stokes, NDIM>::rho()) *=
              1.0 / (1.0 - vars(cellId, LBMVariables<LBEquationType::Navier_Stokes, NDIM>::velocity(dir)));
        }
      }
    }
  } else {
    TERMM(-1, "Invalid LBEquation");
  }
}


#endif // LBM_MOMENTS_H
