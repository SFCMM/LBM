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
/// \param deltaT time
/// \return Velocity at the given time deltaT
template <GInt NDIM, LPTType P>
static constexpr auto freefall_noDrag_vel(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb, const GDouble deltaT)
    -> Point<NDIM> {
  return amb.m_gravity * (1 - amb.m_rho / part.density()) * deltaT + part.initV();
}


/// Analytical solution of the position of a free falling spherical particle without drag.
/// Note: In quiescent flow or no coupling i.e. u_a = 0 (=> No drag)
/// \tparam NDIM Dimensionality
/// \tparam P Particletype
/// \param part Data of a single particle
/// \param amb Ambient condition properties
/// \param deltaT time
/// \return Particle position at the given time deltaT
template <GInt NDIM, LPTType P>
static constexpr auto freefall_noDrag_pos(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb, const GDouble deltaT)
    -> VectorD<NDIM> {
  return 0.5 * amb.m_gravity * (1 - amb.m_rho / part.density()) * gcem::pow(deltaT, 2) + deltaT * part.initV() + part.start();
}


/// Solution to a free falling spherical particle with stokes drag (Re << 1).
/// \tparam NDIM Dimensionality
/// \tparam P Particletype
/// \param part Data of a single particle
/// \param amb Ambient condition properties
/// \param deltaT time
/// \return Velocity at the given time deltaT
template <GInt NDIM, LPTType P>
static constexpr auto freefall_stokes_vel(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb, const GDouble deltaT)
    -> VectorD<NDIM> {
  const GDouble c1 = amb.m_rho / part.density();
  const GDouble c2 = 18.0 * amb.m_nu / (4.0 * part.radius() * part.radius() * part.density()); // Pas /(m^2 * kg/m^3) -> kg/ms / (kg/m) =
                                                                                               // 1/s
  return (gcem::exp(-c2 * deltaT)
          * (c2 * (amb.m_v_infty * (gcem::exp(c2 * deltaT) - 1) + part.initV()) - (c1 - 1) * amb.m_gravity * (gcem::exp(c2 * deltaT) - 1)))
         / c2;
}

/// Solution to a free falling spherical particle with stokes drag (Re << 1).
/// \tparam NDIM Dimensionality
/// \tparam P Particletype
/// \param part Data of a single particle
/// \param amb Ambient condition properties
/// \param deltaT time
/// \return Particle position at the given time deltaT
template <GInt NDIM, LPTType P>
static constexpr auto freefall_stokes_pos(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb, const GDouble deltaT)
    -> VectorD<NDIM> {
  // todo: simplify integration call
  //  auto f_0 = [&](const GDouble _t){
  //    return freefall_stokes_vel(part, amb, _t)[0];
  //  };
  //  auto f_1 = [&](const GDouble _t){
  //    return freefall_stokes_vel(part, amb, _t)[1];
  //  };
  //  auto f_2 = [&](const GDouble _t){
  //    return freefall_stokes_vel(part, amb, _t)[2];
  //  };
  //
  //  // numerical integration
  //  using namespace integration;
  //  multi_simpson(range({0, deltaT}), f_0);
  //  multi_simpson(range({0, deltaT}), f_1);
  //  multi_simpson(range({0, deltaT}), f_2);

  const GInt    num_splits            = 100;
  const GDouble integration_step_size = 0.5 * (deltaT - 0) / static_cast<GDouble>(num_splits);
  VectorD<NDIM> result;
  result.fill(0);

  VectorD<NDIM> f1;
  VectorD<NDIM> f2;
  VectorD<NDIM> f3 = freefall_stokes_vel<NDIM, P>(part, amb, 0);
  for(GInt i = 0; i < num_splits; ++i) {
    f1 = f3;
    f2 = freefall_stokes_vel<NDIM, P>(part, amb, (2 * i + 1) * integration_step_size);
    f3 = freefall_stokes_vel<NDIM, P>(part, amb, (2 * i + 2) * integration_step_size);
    result += 1.0 / 3.0 * integration_step_size * (f1 + 4 * f2 + f3);
  }
  return part.start() + result;
}

/// Solution of the terminal velocity of a spherical particle with stokes drag (Re << 1).
/// \tparam NDIM Dimensionality
/// \tparam P Particletype
/// \param part Data of a single particle
/// \param amb Ambient condition properties
/// \return Terminal velocity of a spherical particle
template <GInt NDIM, LPTType P>
static constexpr auto terminal_velocity_stokes(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb,
                                               const GDouble /*t*/) -> VectorD<NDIM> {
  return amb.m_gravity * gcem::pow(2.0 * part.radius(), 2.0) / (18.0 * amb.m_nu) * part.density();
}

/// Solution of the terminal velocity of a spherical particle with stokes drag (Re << 1) considering buoyancy.
/// \tparam NDIM Dimensionality
/// \tparam P Particletype
/// \param part Data of a single particle
/// \param amb Ambient condition properties
/// \return Terminal velocity of a spherical particle with buoyancy
template <GInt NDIM, LPTType P>
static constexpr auto terminal_velocity_buo_stokes(const ParticleData<NDIM, P>& part, const lpt::AmbientProperties<NDIM>& amb,
                                                   const GDouble /*t*/) -> VectorD<NDIM> {
  return amb.m_gravity * gcem::pow(2.0 * part.radius(), 2.0) / (18.0 * amb.m_nu) * (part.density() - amb.m_rho);
}

} // namespace lpt


/// Obtain the analytical solution by the defined name.
/// \tparam NDIM Dimensionality
/// \tparam P Particletype
/// \param name Name of the analytical solution
/// \return Functional of the analytical solution.
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
