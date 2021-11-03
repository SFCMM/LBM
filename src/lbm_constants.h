#ifndef LBM_LBM_CONSTANTS_H
#define LBM_LBM_CONSTANTS_H

enum class LBInitialConditions {
  MEI06 // as per https://doi.org/10.1016/j.compfluid.2005.08.008
};

enum class LBSolverType {
  BGK,              // classic Bhatnagar–Gross–Krook
  BGK_SMAGORINSKY,  // classic Bhatnagar–Gross–Krook with Smagorinsky SGS model
  BGKC,             // classic Bhatnagar–Gross–Krook-Cercignani (improved accuracy for heat conduction)
  BGKC_SMAGORINSKY, // classic Bhatnagar–Gross–Krook-Cercignani (improved accuracy for heat conduction) with Smagorinsky SGS model
  BGKI,             // classic incompressible Bhatnagar–Gross–Krook
  BGKI_SMAGORINSKY, // classic incompressible Bhatnagar–Gross–Krook with Smagorinsky SGS model
  BGK_THERMAL,      // thermal Bhatnagar–Gross–Krook
  BGK_MRT,          // BGK with multiple relaxation times (improved stability)
  BGK_CDLB,         // cascaded Digital BGK (viscosity is independent of MRT params)
  BGK_CLB,          // cascaded BGK (viscosity is independent of MRT params)
  BGK_CUMULANT,     // cumulant BGK
  BGK_KBC,          // entropy modelled BGK (stable for High Reynolds number on underresolved grids, viscosity is independent of MRT params)
  EULER,
  PE
};


enum class LBMDir { mX, pX, mY, pY, mZ, pZ };

template <GInt NDIM>
inline auto inDirection(const std::array<GDouble, NDIM>& normal, const std::array<GDouble, NDIM>& direction) -> GBool {
  return inDirection(VectorD<NDIM>(&normal[0]), direction);
}

template <GInt NDIM>
inline auto inDirection(const VectorD<NDIM>& normal, const std::array<GDouble, NDIM>& direction) -> GBool {
  return static_cast<GBool>(normal.dot(VectorD<NDIM>(&direction[0])) >= GDoubleEps);
}

enum class BndryType {
  Wall_BounceBack,                   // BounceBack Boundary condition 1st order accurate
  Wall_BounceBack_TangentialVelocity // BounceBack Boundary condition 1st order accurate with tangential velocity
};

enum class LBMethodType { D1Q3, D2Q5, D2Q9, D3Q15, D4Q20, INVALID };

template <GInt NDIM, GInt NDIST>
static constexpr auto getLBMethodType() -> LBMethodType {
  switch(NDIM) {
    case 1:
      return LBMethodType::D1Q3;
    case 2:
      switch(NDIST) {
        case 5:
          return LBMethodType::D2Q5;
        case 9:
          return LBMethodType::D2Q9;
          //        default:
          //          return LBMethodType::INVALID;
      }
    case 3:
      return LBMethodType::D3Q15;
      //    default:
      //      return LBMethodType::INVALID;
    case 4:
      return LBMethodType::D4Q20;
  }
}

static constexpr auto noDists(LBMethodType type) -> GInt {
  switch(type) {
    case LBMethodType::D1Q3:
      return 3;
    case LBMethodType::D2Q5:
      return 5;
    case LBMethodType::D2Q9:
      return 9;
    case LBMethodType::D3Q15:
      return 15;
    case LBMethodType::INVALID:
    default:
      return 0;
  }
}

static constexpr auto dim(LBMethodType type) -> GInt {
  switch(type) {
    case LBMethodType::D1Q3:
      return 1;
    case LBMethodType::D2Q5:
    case LBMethodType::D2Q9:
      return 2;
    case LBMethodType::D3Q15:
      return 3;
    case LBMethodType::INVALID:
    default:
      return 0;
  }
}

template <LBMethodType LBTYPE>
class LBMethod {
 public:
  static constexpr std::array<std::array<GDouble, dim(LBTYPE)>, noDists(LBTYPE)> m_dirs{};
  static constexpr auto oppositeDist(const GInt dist) -> GInt { return 0;}

};


template <>
class LBMethod<LBMethodType::D2Q9> {
 public:
  static constexpr std::array<std::array<GDouble, 2>, 9> m_dirs = {
      {{{-1, 0}}, {{1, 0}}, {{0, -1}}, {{0, 1}}, {{1, 1}}, {{1, -1}}, {{-1, -1}}, {{-1, 1}}, {{0, 0}}}};

  static constexpr auto oppositeDist(const GInt dist) -> GInt {
    switch(dist) {
      case 0:
        return 1;
      case 1:
        return 0;
      case 2:
        return 3;
      case 3:
        return 2;
      case 4:
        return 6;
      case 5:
        return 7;
      case 6:
        return 4;
      case 7:
        return 5;
      case 8:
        return 8;
      default:
        TERMM(-1, "Invalid dist!");
    }
  }
};


#endif // LBM_LBM_CONSTANTS_H
