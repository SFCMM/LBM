#ifndef LBM_ANALYTICAL_SOLUTIONS_H
#define LBM_ANALYTICAL_SOLUTIONS_H

#include <sfcmm_common.h>
#include "lpt/particle.h"
template <GInt NDIM>
using Point = VectorD<NDIM>;

namespace analytical {
struct SolutionConfig {
  GDouble nu = 0;
  GDouble dp = 0;
};

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

// inline auto poiseuille2D(const GDouble ytop, const GDouble ybot, const Point<2> coord, const GDouble maxV) -> Point<2> {
//   return Point<2>({-4.0 * maxV / (gcem::pow(ytop - ybot, 2)) * (coord[1] - ybot) * (coord[1] - ytop), 0});
// }
//
// inline auto poiseuille2D_1(const Point<2> coord) -> Point<2> {
//   //  constexpr GDouble reynoldsNum = 0.7;
//   //  constexpr GDouble relaxTime   = gcem::sqrt(3.0 / 16.0) + 0.5;
//   //  constexpr GDouble refL        = 1.0;
//
//   //  constexpr GDouble dynViscosity = (2.0 * relaxTime - 1.0) / 6.0;
//   // from reynolds obtain reference velocity
//   //  constexpr GDouble refV = reynoldsNum * dynViscosity / refL;
//   constexpr GDouble refV = 0.1;
//   return poiseuille2D(1, 0, coord, refV);
// }

inline auto poiseuille2D(const GDouble ytop, const GDouble ybot, const Point<2> coord, const GDouble dp, const GDouble nu) -> Point<2> {
  return Point<2>({dp / (2.0 * nu) * coord[1] * (ytop - ybot - coord[1]), 0});
}


inline auto poiseuille2D_1(const Point<2> coord, const SolutionConfig* _config = nullptr) -> Point<2> {
  static SolutionConfig config = _config != nullptr ? *_config : config;
  return poiseuille2D(1, 0, coord, config.dp, config.nu);
}

inline auto poiseuille2D_1s(const Point<2> coord) -> Point<2> { return poiseuille2D_1(coord); }


} // namespace ns

namespace euler {}

namespace poisson {
// todo: add reference
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

inline auto poissonSimpleDiffReaction(const Point<1> x) -> Point<1> {
  static constexpr GDouble th = 1.0;
  return Point<1>(gcem::cosh(th * (1.0 - x[0])) / gcem::cosh(th));
}

} // namespace poisson

template <GInt NDIM>
auto getAnalyticalSolution(const GString& name, const SolutionConfig& config) -> std::function<Point<NDIM>(Point<NDIM>)> {
  if constexpr(NDIM == 1) {
    if(name == "poissonCHAI08_1") {
      return &poisson::poissonCHAI08_1;
    }
    if(name == "poissonSimpleDiffReaction") {
      return &poisson::poissonSimpleDiffReaction;
    }
  }

  if constexpr(NDIM == 2) {
    if(name == "couette2D_1_5") {
      return &ns::couette2D_1_5;
    }

    if(name == "poiseuille2D_1") {
      // set config object
      ns::poiseuille2D_1({0, 0}, &config);
      return &ns::poiseuille2D_1s;
    }
  }

  TERMM(-1, "Invalid analyticalSolution :" + name + " selected!");
}

enum class ANALYTICAL_CASE_INDEX { poissonCHAI08_1, poissonSimpleDiffReaction, couette2D, poiseuille2D };
static constexpr std::array<std::string_view, 4> analyticalNames{"poissonCHAI08_1", "poissonSimpleDiffReaction", "couette2D_1_5",
                                                                 "poiseuille2D_1"};

static constexpr auto getStr(ANALYTICAL_CASE_INDEX caseId) -> std::string_view { return analyticalNames[static_cast<GInt>(caseId)]; }

} // namespace analytical

#endif // LBM_ANALYTICAL_SOLUTIONS_H
