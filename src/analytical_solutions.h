#ifndef LBM_ANALYTICAL_SOLUTIONS_H
#define LBM_ANALYTICAL_SOLUTIONS_H

namespace analytical {

namespace ns {
static constexpr auto couette2D(const GDouble couette_wallV, const GDouble couette_channelHeight, const GDouble y) -> GDouble {
  return couette_wallV / couette_channelHeight * y;
}

static constexpr auto couette2D_1_5(const GDouble y) -> GDouble {
  // todo: load properties or compare them???
  constexpr GDouble reynoldsNum = 0.75;
  constexpr GDouble relaxTime   = 0.9;
  constexpr GDouble refL        = 1.0;

  constexpr GDouble dynViscosity = (2.0 * relaxTime - 1.0) / 6.0;
  // from reynolds obtain reference velocity
  constexpr GDouble refV = reynoldsNum * dynViscosity / refL;

  constexpr GDouble couette_top_wallV     = refV;
  constexpr GDouble couette_channelHeight = 5.0;
  return couette2D(couette_top_wallV, couette_channelHeight, y);
}

static constexpr auto poiseuille2D(const GDouble ytop, const GDouble ybot, const GDouble y, const GDouble maxV) -> GDouble {
  return -4 * maxV / (gcem::pow(ytop - ybot, 2)) * (y - ybot) * (y - ytop);
}

static constexpr auto poiseuille2D_1(const GDouble y) -> GDouble {
  constexpr GDouble reynoldsNum = 0.7;
  constexpr GDouble relaxTime   = gcem::sqrt(3.0 / 16.0) + 0.5;
  constexpr GDouble refL        = 1.0;


  constexpr GDouble dynViscosity = (2.0 * relaxTime - 1.0) / 6.0;
  // from reynolds obtain reference velocity
  constexpr GDouble refV = reynoldsNum * dynViscosity / refL;
  return poiseuille2D(1, 0, y, refV);
}
} // namespace ns

namespace euler {}

namespace poisson {}

} // namespace analytical

#endif // LBM_ANALYTICAL_SOLUTIONS_H