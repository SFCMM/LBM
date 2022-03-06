#ifndef LPT_ANALYTICAL_SOLUTIONS_H
#define LPT_ANALYTICAL_SOLUTIONS_H

#include <sfcmm_common.h>
#include "lpt/particle.h"
template <GInt NDIM>
using Point = VectorD<NDIM>;

namespace analytical {
namespace lpt {
template <GInt NDIM>
struct AmbientProperties {
  VectorD<NDIM> m_gravity;
  VectorD<NDIM> m_v_infty;

  GDouble m_rho;
  GDouble m_nu;
};


/// Solution to a free falling spherical particle without drag. In this case the solution is linear!
/// Note: In quiescent flow or no coupling i.e. u_a = 0 (=> No drag)
/// \tparam NDIM Dimensionality
/// \tparam P Particletype
/// \param part Data of a single particle
/// \param amb Ambient condition properties
/// \param t time
/// \return Velocity at the given time t
template <GInt NDIM, LPTType P>
static constexpr auto freefall_noDrag_vel(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb, const GDouble t)
    -> Point<NDIM> {
  return amb.m_gravity * (1 - amb.m_rho / part.density()) * t + part.initV();
}

template <GInt NDIM, LPTType P>
static constexpr auto freefall_noDrag_pos(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb, const GDouble t)
    -> VectorD<NDIM> {
  return 0.5 * amb.m_gravity * (1 - amb.m_rho / part.density()) * gcem::pow(t, 2) + t * part.initV() + part.start();
}

template <GInt NDIM, LPTType P>
static constexpr auto freefall_stokes_vel(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb, const GDouble t)
    -> VectorD<NDIM> {
  const GDouble c1 = amb.m_rho / part.density();
  const GDouble c2 = 18.0 * amb.m_nu / (4.0 * part.radius() * part.radius() * part.density()); // Pas /(m^2 * kg/m^3) -> kg/ms / (kg/m) =
                                                                                               // 1/s
  return (gcem::exp(-c2 * t)
          * (c2 * (amb.m_v_infty * (gcem::exp(c2 * t) - 1) + part.initV()) - (c1 - 1) * amb.m_gravity * (gcem::exp(c2 * t) - 1)))
         / c2;
}

template <GInt NDIM, LPTType P>
static constexpr auto freefall_stokes_pos(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb, const GDouble t)
    -> VectorD<NDIM> {
  const GDouble c1 = amb.m_rho / part.density();
  const GDouble c2 =
      18.0 * amb.m_nu / (4.0 * part.radius() * part.radius() * part.density()); // Pas /(m^2 * kg/m^3) -> kg/ms / (kg/m) = 1/s
  return (gcem::exp(-c2 * t) * (-c1 * amb.m_gravity + c2 * amb.m_v_infty - c2 * part.initV() + amb.m_gravity)) / (c2 * c2)
         + (t * (-c1 * amb.m_gravity + c2 * amb.m_v_infty + amb.m_gravity)) / c2 + part.start();
}

template <GInt NDIM, LPTType P>
static constexpr auto terminal_velocity_stokes(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb,
                                               const GDouble /*t*/) -> VectorD<NDIM> {
  return amb.m_gravity * gcem::pow(2.0 * part.radius(), 2.0) / (18.0 * amb.m_nu) * part.density();
}

template <GInt NDIM, LPTType P>
static constexpr auto terminal_velocity_buo_stokes(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb,
                                                   const GDouble /*t*/) -> VectorD<NDIM> {
  return amb.m_gravity * gcem::pow(2.0 * part.radius(), 2.0) / (18.0 * amb.m_nu) * (part.density() - amb.m_rho);
}

// template <GInt NDIM, LPTType P>
// static constexpr auto terminal_velocity_buo_stokes(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb,
//                                                    const GDouble /*t*/) -> VectorD<NDIM> {
//   const GDouble Re = amb.m_rho * 2 * part.radius() * part.velocity().norm() / amb.m_nu;
//   // stokes drag
//   const GDouble CD = drag();
//   if(amb.m_rho < part.density()) {
//     return gcem::sqrt(4 * amb.m_gravity / (3 * CD) * (part.density() - amb.m_rho) / amb.m_rho);
//   }
//   // particle has a lower density than the surrounding medium so it moves upwards (v<0)
//   return -gcem::sqrt(4 * amb.m_gravity / (3 * -CD) * (part.density() - amb.m_rho) / amb.m_rho);
// }

} // namespace lpt

template <GInt NDIM, LPTType P>
auto getAnalyticalSolution(const GString& name)
    -> std::function<VectorD<NDIM>(const ParticleData<NDIM, P>&, const lpt::AmbientProperties<NDIM>&, const GDouble)> {
  if(name == "freefall_stokes_vel") {
    return &lpt::freefall_stokes_vel<NDIM, P>;
  }
  if(name == "freefall_stokes_pos") {
    return &lpt::freefall_stokes_pos<NDIM, P>;
  }

  if(name == "freefall_nodrag_vel") {
    return &lpt::freefall_noDrag_vel<NDIM, P>;
  }
  if(name == "freefall_nodrag_pos") {
    return &lpt::freefall_noDrag_pos<NDIM, P>;
  }

  if(name == "terminal_stokes_vel") {
    return &lpt::terminal_velocity_stokes<NDIM, P>;
  }
  if(name == "terminal_stokes_buo_vel") {
    return &lpt::terminal_velocity_buo_stokes<NDIM, P>;
  }


  TERMM(-1, "Invalid analyticalSolution selected! Selected type: " + name);
}

} // namespace analytical

#endif // LPT_ANALYTICAL_SOLUTIONS_H
