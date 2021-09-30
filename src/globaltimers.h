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

    GridGeneratorTotal,
    GridGeneration,
    GridInit,
    GridPart,
    GridUniform,
    GridRefinement,
    GridIo,

    LBMSolverTotal,

    // counter
    _count
  };
};
} // namespace gridgenerator

inline std::array<GInt, gridgenerator::Timers_::_count> TimeKeeper{}; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
using Timers = gridgenerator::Timers_;

#endif // GRIDGENERATOR_GLOBALTIMERS_H
