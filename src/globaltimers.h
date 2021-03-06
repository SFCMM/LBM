#ifndef GRIDGENERATOR_GLOBALTIMERS_H
#define GRIDGENERATOR_GLOBALTIMERS_H

namespace gridgenerator {
// Create struct for easy timer identification
struct Timers_ {
  // Enum to store timer "names"
  enum {
    // Timer groups
    AppGroup,


    timertotal,
    IO,
    Init,

    // Grid Generator
    GridGeneratorTotal,
    GridGeneration,
    GridInit,
    GridPart,
    GridUniform,
    GridRefinement,
    GridIo,

    // LBM
    LBMSolverTotal,
    LBMInit,
    LBMMainLoop,
    LBMCalc,
    LBMProp,
    LBMColl,
    LBMBnd,
    LBMEq,
    LBMMacro,
    LBMForce,
    LBMPost,
    LBMIo,

    // Particle
    LPTSolverTotal,
    LPTInit,
    LPTMainLoop,
    LPTCalc,
    LPTInt,
    LPTColl,
    LPTForce,
    LPTPost,
    LPTIo,

    // counter
    _count
  };
};
} // namespace gridgenerator

inline std::array<GInt, gridgenerator::Timers_::_count> TimeKeeper{}; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
using Timers = gridgenerator::Timers_;

#endif // GRIDGENERATOR_GLOBALTIMERS_H
