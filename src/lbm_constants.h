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

enum class LBMethodType { D1Q3, D2Q5, D2Q9 };

namespace LBMethod {
static constexpr std::array<GInt, 3> noDists = {3, 5, 9};
static constexpr std::array<GInt, 3> dim     = {3, 5, 9};
} // namespace LBMethod

#endif // LBM_LBM_CONSTANTS_H
