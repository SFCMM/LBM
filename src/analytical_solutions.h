#ifndef LBM_ANALYTICAL_SOLUTIONS_H
#define LBM_ANALYTICAL_SOLUTIONS_H

namespace analytical {

namespace ns {
static constexpr auto couette2D(const GDouble couette_wallV, const GDouble couette_channelHeight, const GDouble y) -> GDouble {
  return couette_wallV / couette_channelHeight * y;
}

static constexpr auto couette2D_1_5(const GDouble y) {
  // todo: load properties or compare them???
  constexpr GDouble reynoldsNum  = 0.75;
  constexpr GDouble relaxTime    = 0.9;
  constexpr GDouble dynViscosity = (2.0 * relaxTime - 1.0) / 6.0;
  constexpr GDouble refL         = 1.0;

  // from reynolds obtain reference velocity
  constexpr GDouble refV = reynoldsNum * dynViscosity / refL;

  constexpr GDouble couette_top_wallV     = refV;
  constexpr GDouble couette_channelHeight = 5.0;
  return couette2D(couette_top_wallV, couette_channelHeight, y);
}
} // namespace ns

namespace euler {}

namespace poisson {}

} // namespace analytical

#endif // LBM_ANALYTICAL_SOLUTIONS_H
