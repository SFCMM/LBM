#ifndef LBM_PV_H
#define LBM_PV_H
namespace PV {
static constexpr auto U() -> GInt { return 0; }

static constexpr auto V() -> GInt { return 1; }

static constexpr auto W() -> GInt { return 2; }

template <GInt NDIM>
static constexpr auto velocities() -> std::array<GInt, NDIM> {
  if constexpr(NDIM == 1) {
    return {U()};
  }
  if constexpr(NDIM == 2) {
    return {U(), V()};
  }
  if constexpr(NDIM == 3) {
    return {U(), V(), W()};
  }
  TERMM(-1, "Not implemented!");
}

template <GInt NDIM>
static constexpr auto rho() -> GInt {
  return velocities<NDIM>()[NDIM - 1] + 1;
}

template <GInt NDIM>
static constexpr auto temperature() -> GInt {
  return rho<NDIM>() + 1;
}

} // namespace PV


#endif // LBM_PV_H
