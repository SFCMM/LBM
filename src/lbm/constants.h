#ifndef LBM_CONSTANTS_H
#define LBM_CONSTANTS_H

enum class LBMethodType;
template <LBMethodType LBTYPE>
class LBMethod;

// default environment properties
static constexpr GDouble defaultT20C       = 293.15;
static constexpr GDouble defaultMachNumber = 0.01;

// default numerical properties
static constexpr GDouble defaultRelaxT = 0.9;

// default solver setting
static constexpr GInt defaultInfoOutInterval  = 10;
static constexpr GInt defaultSolutionInterval = 100;
static constexpr GInt maxNumberDistributions  = 30;

// default speed of sound
static constexpr GDouble lbm_cs   = 1.0 / gcem::sqrt(3.0); // assuming dx/dt = 1
static constexpr GDouble lbm_cssq = 1.0 / 3.0;             // assuming dx/dt = 1

enum class LBEquation { Navier_Stokes, Poisson, Navier_Stokes_Poisson };
static constexpr std::array<std::string_view, 3> LBEquationName = {"Navier-Stokes",
                                                                   // Implementation based on "A novel lattice Boltzmann model for the
                                                                   // Poisson equation", 2008, Applied Mathematical Modeling, Chai and
                                                                   // Shi [CHAI08]
                                                                   "Poisson", "Navier-Stokes-Poisson"};

static constexpr auto getLBEquationType(const std::string_view equationName) -> LBEquation {
  if(equationName == "navierstokes") {
    return LBEquation::Navier_Stokes;
  }
  if(equationName == "poisson") {
    return LBEquation::Poisson;
  }
  if(equationName == "navierstokespoisson") {
    return LBEquation::Navier_Stokes_Poisson;
  }
  TERMM(-1, "Invalid equation configuration!");
}

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

template <GInt NDIM>
inline auto inDirection(const std::array<GDouble, NDIM>& normal, const std::array<GDouble, NDIM>& direction) -> GBool {
  return inDirection(VectorD<NDIM>(normal.data()), direction);
}

template <GInt NDIM>
inline auto inDirection(const VectorD<NDIM>& normal, const std::array<GDouble, NDIM>& direction) -> GBool {
  return static_cast<GBool>(normal.dot(VectorD<NDIM>(direction.data())) >= GDoubleEps);
}

enum class BndryType {
  Wall_BounceBack,                    // BounceBack Wall Boundary condition 2nd order accurate but viscosity depends on relaxation
  Wall_BounceBack_TangentialVelocity, // BounceBack Boundary condition 2nd order accurate with tangential velocity but viscosity depends
                                      // on relaxation
  Wall_Equilibrium,                   // Wet Node Equilibrium method 1st order accurate
  Wall_NEEM,                          // Wet Node Non-equilibrium extrapolation method 2nd order accurate
  Wall_NEBB,                          // Wet Node Non-equilibrium method bounceback method 3rd order accurate
  Inlet_BounceBack_ConstPressure,     // BounceBack Boundary condition 1st order accurate with tangential velocity
  Outlet_BounceBack_ConstPressure,    // BounceBack Boundary condition 1st order accurate with tangential velocity
  Periodic,                           // Periodic boundary condition handled through boundary class
  Dirichlet_NEEM                      // Dirichlet boundary condition using Non-equilibrium extrapolation method
};

enum class LBMethodType { D1Q3, D2Q5, D2Q9, D3Q19, D3Q27, D4Q40, INVALID };


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
        default:
          std::terminate();
      }
    case 3:
      switch(NDIST) {
        case 19:
          return LBMethodType::D3Q19;
        case 27:
          return LBMethodType::D3Q27;
        default:
          std::terminate();
      }
    case 4:
      return LBMethodType::D4Q40;
    default:
      std::terminate();
  }
}

static constexpr auto getLBMethodType(const std::string_view modelName) -> LBMethodType {
  if(modelName == "D1Q3") {
    return LBMethodType::D1Q3;
  }
  if(modelName == "D2Q5") {
    return LBMethodType::D2Q5;
  }
  if(modelName == "D2Q9") {
    return LBMethodType::D2Q9;
  }
  if(modelName == "D3Q19") {
    return LBMethodType::D3Q19;
  }
  if(modelName == "D3Q27") {
    return LBMethodType::D3Q27;
  }
  TERMM(-1, "Invalid model configuration!");
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
        default:
          std::terminate();
      }
    case 3:
      switch(noDistributions) {
        case 19:
          return getLBMethodType<3, 19>();
        case 27:
          return getLBMethodType<3, 27>();
        default:
          std::terminate();
      }
    case 4:
      return getLBMethodType<4, 20>();
    default:
      std::terminate();
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
    case LBMethodType::D3Q19:
      return 19;
    case LBMethodType::D3Q27:
      return 27;
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
    case LBMethodType::D3Q19:
    case LBMethodType::D3Q27:
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
  static constexpr auto                                                          oppositeDist(const GInt /*dist*/) -> GInt {
    TERMM(-1, "Invalid call");
    return -1;
  }
  static constexpr GInt                   m_dim     = 0;
  static constexpr GInt                   m_noDists = 0;
  static constexpr std::array<GDouble, 1> m_weights = {0};

  static constexpr GDouble                m_poissonAlpha   = NAN;
  static constexpr std::array<GDouble, 1> m_poissonWeights = {0};

  static constexpr GBool            m_isThermal  = false;
  static constexpr GBool            m_canPoisson = false;
  static constexpr std::string_view m_name       = "INVALID";
};

template <>
class LBMethod<LBMethodType::D1Q3> {
 public:
  static constexpr std::array<std::array<GDouble, 1>, 3> m_dirs = {{{{-1}}, {{1}}, {{0}}}};

  static constexpr std::array<GInt, 3> m_oppositeDist = {1, 0, 2};

  /// look-up table for opposite direction
  static constexpr auto oppositeDist(const GInt dist) -> GInt { return m_oppositeDist[dist]; }

  static constexpr std::array<GDouble, 3> m_weights = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};

  // see equation 2.3 in CHAI08
  static constexpr GDouble m_poissonAlpha = 1.0 / 3.0;
  // see CHAI08
  static constexpr std::array<GDouble, 3> m_poissonWeights = {0.5, 0.5, 0};

  static constexpr GInt             m_dim        = 1;
  static constexpr GInt             m_noDists    = 3;
  static constexpr GBool            m_isThermal  = false;
  static constexpr GBool            m_canPoisson = true;
  static constexpr std::string_view m_name       = "D1Q3";
};

template <>
class LBMethod<LBMethodType::D2Q5> {
 public:
  static constexpr std::array<std::array<GDouble, 2>, 5> m_dirs = {{{{-1, 0}}, {{1, 0}}, {{0, -1}}, {{0, 1}}, {{0, 0}}}};

  static constexpr std::array<GInt, 5> m_oppositeDist = {1, 0, 3, 2, 4};

  /// look-up table for opposite direction
  static constexpr auto oppositeDist(const GInt dist) -> GInt { return m_oppositeDist[dist]; }

  static constexpr std::array<GDouble, 5> m_weights = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0};

  // see equation 2.3 in CHAI08
  static constexpr GDouble m_poissonAlpha = 2.0 / 3.0; // test gave 1/3
  // see CHAI08
  static constexpr std::array<GDouble, 5> m_poissonWeights = {0.25, 0.25, 0.25, 0.25, 0.0};

  static constexpr GInt             m_dim        = 2;
  static constexpr GInt             m_noDists    = 5;
  static constexpr GBool            m_isThermal  = false;
  static constexpr GBool            m_canPoisson = true;
  static constexpr std::string_view m_name       = "D2Q5";
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

  static constexpr GDouble                m_poissonAlpha   = 1.0 / 3.0;
  static constexpr std::array<GDouble, 9> m_poissonWeights = {1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0,
                                                              1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 0.0};


  static constexpr GInt             m_dim        = 2;
  static constexpr GInt             m_noDists    = 9;
  static constexpr GBool            m_isThermal  = false;
  static constexpr GBool            m_canPoisson = false;
  static constexpr std::string_view m_name       = "D2Q9";
};

template <>
class LBMethod<LBMethodType::D3Q19> {
 public:
  static constexpr std::array<std::array<GDouble, 3>, 19> m_dirs = {
      {{{-1, 0, 0}},
       {{1, 0, 0}}, // main dir x
       {{0, -1, 0}},
       {{0, 1, 0}}, // main dir y
       {{0, 0, -1}},
       {{0, 0, 1}}, // main dir z
       {{-1, -1, 0}},
       {{-1, 1, 0}},
       {{1, -1, 0}},
       {{1, 1, 0}}, // diagonals in 2D plane (Note: order is different to 2DQ9)
       {{-1, 0, -1}},
       {{-1, 0, 1}},
       {{1, 0, -1}},
       {{1, 0, 1}}, // diagonals in x dir
       {{0, -1, -1}},
       {{0, -1, 1}},
       {{0, 1, -1}},
       {{0, 1, 1}},   // diagonals in y dir
       {{0, 0, 0}}}}; // center

  static constexpr std::array<GInt, 19> m_oppositeDist = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 18};

  /// look-up table for opposite direction
  static constexpr auto oppositeDist(const GInt dist) -> GInt { return m_oppositeDist[dist]; }

  static constexpr std::array<GDouble, 19> m_weights = {1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 36.0,
                                                        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
                                                        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 3.0};

  static constexpr GDouble m_poissonAlpha = NAN;
  //  static constexpr std::array<GDouble, 19> m_poissonWeights = {1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0,
  //                                                              1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 0.0};


  static constexpr GInt             m_dim        = 3;
  static constexpr GInt             m_noDists    = 19;
  static constexpr GBool            m_isThermal  = false;
  static constexpr GBool            m_canPoisson = false;
  static constexpr std::string_view m_name       = "D3Q19";
};

template <>
class LBMethod<LBMethodType::D3Q27> {
 public:
  static constexpr std::array<std::array<GDouble, 3>, 27> m_dirs = {{{{-1, 0, 0}},
                                                                     {{1, 0, 0}}, // main dir x
                                                                     {{0, -1, 0}},
                                                                     {{0, 1, 0}}, // main dir y
                                                                     {{0, 0, -1}},
                                                                     {{0, 0, 1}}, // main dir z
                                                                     {{-1, -1, 0}},
                                                                     {{-1, 1, 0}},
                                                                     {{1, -1, 0}},
                                                                     {{1, 1, 0}}, // diagonals in 2D plane (Note: order is different to
                                                                                  // 2DQ9)
                                                                     {{-1, 0, -1}},
                                                                     {{-1, 0, 1}},
                                                                     {{1, 0, -1}},
                                                                     {{1, 0, 1}}, // diagonals in x dir
                                                                     {{0, -1, -1}},
                                                                     {{0, -1, 1}},
                                                                     {{0, 1, -1}},
                                                                     {{0, 1, 1}}, // diagonals in y dir
                                                                     {{-1.0, -1.0, -1.0}},
                                                                     {{-1.0, -1.0, 1.0}},
                                                                     {{-1.0, 1.0, -1.0}},
                                                                     {{-1.0, 1.0, 1.0}},
                                                                     {{1.0, -1.0, -1.0}},
                                                                     {{1.0, -1.0, 1.0}},
                                                                     {{1.0, 1.0, -1.0}},
                                                                     {{1.0, 1.0, 1.0}}, // tridiagonal directions
                                                                     {{0, 0, 0}}}};     // center

  static constexpr std::array<GInt, 27> m_oppositeDist = {1,  0,  3,  2,  5,  4,  9,  8,  7,  6,  13, 12, 11, 10,
                                                          17, 16, 15, 14, 25, 24, 23, 22, 21, 20, 19, 18, 26};

  /// look-up table for opposite direction
  static constexpr auto oppositeDist(const GInt dist) -> GInt { return m_oppositeDist[dist]; }

  static constexpr std::array<GDouble, 27> m_weights = {
      2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
      1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
      1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 8.0 / 27.0};

  static constexpr GDouble m_poissonAlpha = NAN;
  //  static constexpr std::array<GDouble, 19> m_poissonWeights = {1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0,
  //                                                              1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 0.0};


  static constexpr GInt             m_dim        = 3;
  static constexpr GInt             m_noDists    = 27;
  static constexpr GBool            m_isThermal  = false;
  static constexpr GBool            m_canPoisson = false;
  static constexpr std::string_view m_name       = "D3Q27";
};

template <LBMethodType LBTYPE>
static constexpr auto noVars(const LBEquation equation) -> GInt {
  switch(equation) {
    case LBEquation::Navier_Stokes:
      // Velocities + Density + Temperature
      return LBMethod<LBTYPE>::m_dim + 1 + static_cast<GInt>(LBMethod<LBTYPE>::m_isThermal);
    case LBEquation::Poisson:
      // Potential
      return 1;
    case LBEquation::Navier_Stokes_Poisson:
      // Velocities + Density + Temperature + Potential
      return LBMethod<LBTYPE>::m_dim + 1 + static_cast<GInt>(LBMethod<LBTYPE>::m_isThermal) + 1;
    default:
      TERMM(-1, "Invalid equation type!");
  }
}


#endif // LBM_CONSTANTS_H
