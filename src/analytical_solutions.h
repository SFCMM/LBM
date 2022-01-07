#ifndef LBM_ANALYTICAL_SOLUTIONS_H
#define LBM_ANALYTICAL_SOLUTIONS_H

#include <sfcmm_common.h>
#include "lpt/particle.h"
template <GInt NDIM>
using Point = VectorD<NDIM>;

namespace analytical {

namespace ns {
inline auto couette2D(const GDouble couette_wallV, const GDouble couette_channelHeight, const Point<2> coord) -> Point<2> {
  return Point<2>({couette_wallV / couette_channelHeight * coord[1], 0});
}

inline auto couette2D_1_5(const Point<2> coord) -> Point<2> {
  // todo: load properties or compare them???
  constexpr GDouble reynoldsNum = 0.75;
  constexpr GDouble relaxTime   = 0.9;
  constexpr GDouble refL        = 1.0;

  constexpr GDouble dynViscosity = (2.0 * relaxTime - 1.0) / 6.0;
  // from reynolds obtain reference velocity
  constexpr GDouble refV = reynoldsNum * dynViscosity / refL;

  constexpr GDouble couette_top_wallV     = refV;
  constexpr GDouble couette_channelHeight = 5.0;
  return couette2D(couette_top_wallV, couette_channelHeight, coord);
}

inline auto poiseuille2D(const GDouble ytop, const GDouble ybot, const Point<2> coord, const GDouble maxV) -> Point<2> {
  return Point<2>({-4.0 * maxV / (gcem::pow(ytop - ybot, 2)) * (coord[1] - ybot) * (coord[1] - ytop), 0});
}

inline auto poiseuille2D_1(const Point<2> coord) -> Point<2> {
  //  constexpr GDouble reynoldsNum = 0.7;
  //  constexpr GDouble relaxTime   = gcem::sqrt(3.0 / 16.0) + 0.5;
  //  constexpr GDouble refL        = 1.0;

  //  constexpr GDouble dynViscosity = (2.0 * relaxTime - 1.0) / 6.0;
  // from reynolds obtain reference velocity
  //  constexpr GDouble refV = reynoldsNum * dynViscosity / refL;
  constexpr GDouble refV = 0.1;
  return poiseuille2D(1, 0, coord, refV);
}


} // namespace ns

namespace euler {}

namespace poisson {
///
/// \param x Position
/// \return Solution as given by CHAI08
inline auto poissonCHAI08_1(const Point<1> x) -> Point<1> {
  // see equation 3.4 in CHAI08
  constexpr GDouble k       = 27.79; // Debye-Hueckel approximation
  constexpr GDouble exp_pK  = gcem::exp(k);
  constexpr GDouble exp_mK  = gcem::exp(-k);
  GDouble           exp_pKx = gcem::exp(k * x[0]);
  GDouble           exp_mKx = gcem::exp(-k * x[0]);

  return Point<1>((exp_pK - 1.0) / (exp_pK - exp_mK) * exp_mKx + (1.0 - exp_mK) / (exp_pK - exp_mK) * exp_pKx);
}
} // namespace poisson

namespace lpt {
template <GInt NDIM, LPTType P>
static constexpr auto freefall_noDrag_vel(const ParticleData<NDIM, P>& part, const VectorD<NDIM>& init_v, const GDouble t,
                                          const Point<NDIM>& gravity, const GDouble rho_a = 0) -> Point<NDIM> {
  return gravity * (1 - rho_a / part.density()) * t + init_v;
}

template <GInt NDIM, LPTType P>
static constexpr auto freefall_noDrag_pos(const Point<NDIM>& start, const ParticleData<NDIM, P>& part, const VectorD<NDIM>& init_v,
                                          const GDouble t, const Point<NDIM>& gravity, const GDouble rho_a = 0) -> VectorD<NDIM> {
  return 0.5 * gravity * (1 - rho_a / part.density()) * gcem::pow(t, 2) + t * init_v + start;
}

template <GInt NDIM, LPTType P>
static constexpr auto freefall_stokes_vel(const ParticleData<NDIM, P>& part, const VectorD<NDIM>& init_v, const GDouble nu_a,
                                          const GDouble rho_a, const VectorD<NDIM>& v_a, const GDouble t, const VectorD<NDIM>& gravity)
    -> VectorD<NDIM> {
  const GDouble c1 = rho_a / part.density();
  const GDouble c2 = 18.0 * nu_a / (4.0 * part.radius() * part.radius() * part.density()); // Pas /(m^2 * kg/m^3) -> kg/ms / (kg/m) = 1/s
  return (gcem::exp(-c2 * t) * (c2 * (v_a * (gcem::exp(c2 * t) - 1) + init_v) - (c1 - 1) * gravity * (gcem::exp(c2 * t) - 1))) / c2;
}

template <GInt NDIM, LPTType P>
static constexpr auto freefall_stokes_pos(const Point<NDIM>& start, const ParticleData<NDIM, P>& part, const VectorD<NDIM>& init_v,
                                          const GDouble nu_a, const GDouble rho_a, const VectorD<NDIM>& v_a, const GDouble t,
                                          const VectorD<NDIM>& gravity) -> VectorD<NDIM> {
  const GDouble c1 = rho_a / part.density();
  const GDouble c2 = 18.0 * nu_a / (4.0 * part.radius() * part.radius() * part.density()); // Pas /(m^2 * kg/m^3) -> kg/ms / (kg/m) = 1/s
  return (gcem::exp(-c2 * t)(-c1 * gravity + c2 * v_a - c2 * init_v + gravity)) / (c2 * c2)
         + (t * (-c1 * gravity + c2 * v_a + gravity)) / c2 + start;
}
} // namespace lpt

template <GInt NDIM>
auto getAnalyticalSolution(const GString& name) -> std::function<Point<NDIM>(Point<NDIM>)> {
  if constexpr(NDIM == 1) {
    if(name == "poissonCHAI08_1") {
      return &poisson::poissonCHAI08_1;
    }
  }

  if constexpr(NDIM == 2) {
    if(name == "couette2D_1_5") {
      return &ns::couette2D_1_5;
    }

    if(name == "poiseuille2D_1") {
      return &ns::poiseuille2D_1;
    }
  }

  TERMM(-1, "Invalid analyticalSolution selected!");
}

enum class ANALYTICAL_CASE_INDEX { poissonCHAI08_1, couette2D, poiseuille2D };
static constexpr std::array<std::string_view, 3> analyticalNames{"poissonCHAI08_1", "couette2D_1_5", "poiseuille2D_1"};

static constexpr auto getStr(ANALYTICAL_CASE_INDEX caseId) -> std::string_view { return analyticalNames[static_cast<GInt>(caseId)]; }

} // namespace analytical

#endif // LBM_ANALYTICAL_SOLUTIONS_H
