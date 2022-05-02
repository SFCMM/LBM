#ifndef LBM_EQUILIBRIUM_FUNC_H
#define LBM_EQUILIBRIUM_FUNC_H

#include <sfcmm_common.h>
#include "constants.h"


namespace eq {


/// linearised Conservation equation (e.g. for acoustics)
/// \param weight LBM lattice weight
/// \param density density
/// \param density0 density at rest
/// \param ci_alpha
/// \param u_alpha
/// \return equilibrium distribution of the linearised conservation equations
static inline auto acoustic(const GDouble weight, const GDouble density, const GDouble density0, const GDouble ci_alpha,
                            const GDouble u_alpha) -> GDouble {
  // todo: add test
  return weight * (density + density0 * ci_alpha * u_alpha / lbm_cssq);
}

/// Navier-Stokes equation for incompressible *steady-state* flows with Ma << 1
/// \param weight LBM lattice weight
/// \param density Density
/// \param density0 Density at rest
/// \param ci_alpha
/// \param ci_beta
/// \param u_alpha
/// \param u_beta
/// \param delta_alpha_beta
/// \return equilibrium distribution of the incompressible Navier-Stokes
static inline auto incompressible(const GDouble weight, const GDouble density, const GDouble density0, const GDouble ci_alpha,
                                  const GDouble ci_beta, const GDouble u_alpha, const GDouble u_beta, const GDouble delta_alpha_beta)
    -> GDouble {
  // todo: add test
  return weight
         * (density
            + density0
                  * (ci_alpha * u_alpha / lbm_cssq
                     + u_alpha * u_beta * (ci_alpha * ci_beta - lbm_cssq * delta_alpha_beta) / (2.0 * lbm_cssq * lbm_cssq)));
}


/// Normal LBM equilibrium distribution for Navier-Stokes
/// \param weight LBM lattice weight
/// \param density density
/// \param cu c_[i][dir] * u[dir]
/// \param vsq velocity squared (u^2)
/// \return value of the equilibrium distribution function for Navier-Stokes
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

template <LBMethodType LBTYPE>
static inline auto defaultEq(const GInt dist, const GDouble density, const GDouble* const velocity) -> GDouble {
  //  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;

  GDouble vsq = 0;
  for(GInt dir = 0; dir < NDIM; ++dir) {
    vsq += velocity[dir] * velocity[dir];
  }

  GDouble cu = 0;
  for(GInt dir = 0; dir < NDIM; ++dir) {
    cu += velocity[dir] * LBMethod<LBTYPE>::m_dirs[dist][dir];
  }
  return eq::defaultEq(LBMethod<LBTYPE>::m_weights[dist], density, cu, vsq);
}

/// Symmetric LBM equilibrium distribution for Navier-Stokes (used in bnds)
/// \param weight LBM lattice weight
/// \param density density
/// \param cu c_[i][dir] * u[dir]
/// \param vsq velocity squared (u^2)
/// \return value of the equilibrium distribution function for Navier-Stokes
static inline auto symmEq(const GDouble weight, const GDouble density, const GDouble cu, const GDouble vsq) -> GDouble {
  return weight * density * (1.0 + cu * cu / (2.0 * lbm_cssq * lbm_cssq) - vsq / (2.0 * lbm_cssq));
}

template <LBMethodType LBTYPE>
static inline void symmEq(GDouble* feq, const GDouble density, const GDouble* const velocity) {
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
    feq[dist] = eq::symmEq(LBMethod<LBTYPE>::m_weights[dist], density, cu, vsq);
  }
}

template <LBMethodType LBTYPE>
static inline auto symmEq(const GInt dist, const GDouble density, const GDouble* const velocity) -> GDouble {
  //  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;

  GDouble vsq = 0;
  for(GInt dir = 0; dir < NDIM; ++dir) {
    vsq += velocity[dir] * velocity[dir];
  }

  GDouble cu = 0;
  for(GInt dir = 0; dir < NDIM; ++dir) {
    cu += velocity[dir] * LBMethod<LBTYPE>::m_dirs[dist][dir];
  }
  return eq::symmEq(LBMethod<LBTYPE>::m_weights[dist], density, cu, vsq);
}

/// Equilibrium distribution for the poisson equation
/// \param weight LBM Lattice weight
/// \param potential Potential
/// \return Value of the equilibrium function
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

#endif // LBM_EQUILIBRIUM_FUNC_H
