#ifndef LBM_LBM_CONSTANTS_H
#define LBM_LBM_CONSTANTS_H

enum class LBMethodType;
template <LBMethodType LBTYPE>
class LBMethod;

// default environment properties
static constexpr GDouble defaultT20C       = 293.15;
static constexpr GDouble defaultMachNumber = 0.2;

// default numerical properties
static constexpr GDouble defaultRelaxT = 0.9;

// default solver setting
static constexpr GInt defaultInfoOutInterval  = 10;
static constexpr GInt defaultSolutionInterval = 100;
static constexpr GInt maxNumberDistributions  = 30;

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

static constexpr auto getLBMethodType(const GInt noDims, const GInt noDistributions) -> LBMethodType {
  switch(noDims) {
    case 1:
      return getLBMethodType<1, 3>();
    case 2:
      switch(noDistributions) {
        case 5:
          return getLBMethodType<2, 5>();
        case 9:
          return getLBMethodType<2, 9>();
          //        default:
          //          return LBMethodType::INVALID;
      }
    case 3:
      return getLBMethodType<3, 15>();
      //    default:
      //      return LBMethodType::INVALID;
    case 4:
      return getLBMethodType<4, 20>();
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
  static constexpr auto                                                          oppositeDist(const GInt dist) -> GInt { return 0; }
  static constexpr GInt                                                          m_dim     = 0;
  static constexpr GInt                                                          m_noDists = 0;
  static constexpr std::array<GDouble, 1> m_weights = {0};
  static constexpr GBool m_isThermal = false;
};


template <>
class LBMethod<LBMethodType::D2Q9> {
 public:
  static constexpr std::array<std::array<GDouble, 2>, 9> m_dirs = {
      {{{-1, 0}}, {{1, 0}}, {{0, -1}}, {{0, 1}}, {{1, 1}}, {{1, -1}}, {{-1, -1}}, {{-1, 1}}, {{0, 0}}}};

  static constexpr std::array<GInt, 9> m_oppositeDist = {1, 0, 3, 2, 6, 7, 4, 5, 8};

  /// look-up table for opposite direction
  static constexpr auto oppositeDist(const GInt dist) -> GInt { return m_oppositeDist[dist]; }

  static constexpr std::array<GDouble, 9> m_weights = {1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0, 1.0 / 36.0,
                                                       1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 4.0 / 9.0};

  static constexpr GInt m_dim     = 2;
  static constexpr GInt m_noDists = 9;
  static constexpr GBool m_isThermal = false;
};


#endif // LBM_LBM_CONSTANTS_H
