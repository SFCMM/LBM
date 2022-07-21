// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_CONSTANTS_H
#define SFCMM_CONSTANTS_H

#include <array>
#include <limits>
#include <vector>
#include <iostream>
#include "sfcmm_types.h"

static constexpr GInt    BASE2          = 2;
static constexpr GDouble HALF           = 0.5;
static constexpr GInt    INVALID_CELLID = -1;
static constexpr GInt    INVALID_DIR    = -1;

static constexpr GDouble PI = 3.141592653589793238462643383279;

/// Generate invalid list for init of arrays
/// \tparam LENGTH Length of the invalid list
/// \return Invalid list of LENGTH
template <GInt LENGTH>
static constexpr auto INVALID_LIST() -> std::array<GInt, LENGTH> {
  std::array<GInt, LENGTH> invalid{};
  std::fill_n(invalid.begin(), LENGTH, INVALID_CELLID);
  return invalid;
}

/// Generate nan list for init of arrays
/// \tparam LENGTH Length of the nan list
/// \return NAN list of LENGTH
template <GInt LENGTH>
static constexpr auto NAN_LIST() -> std::array<GDouble, LENGTH> {
  std::array<GDouble, LENGTH> invalid{};
  std::fill_n(invalid.begin(), LENGTH, NAN);
  return invalid;
}

// Memory
static constexpr GInt    KBIT  = 1024;
static constexpr GDouble DKBIT = static_cast<GDouble>(KBIT);

// Time
namespace timeconst {
static constexpr GInt    MINUTE  = 60;
static constexpr GDouble DMINUTE = 60;
static constexpr GInt    HOUR    = 3600;
static constexpr GDouble DHOUR   = 3600;
static constexpr GInt    DAY     = HOUR * 24;
static constexpr GDouble DDAY    = DHOUR * 24;
static constexpr GInt    WEEK    = DAY * 7;
static constexpr GDouble DWEEK   = DDAY * 7;
} // namespace timeconst

enum class Debug_Level { no_debug, debug, max_debug };
static constexpr std::array<std::string_view, 3> DEBUG_LEVEL = {"NO DEBUG", "DEBUG", "MAXIMUM DEBUG"};

enum class SolverType { NONE, GRIDDER, LBM, LPT };
static constexpr std::array<std::string_view, 4> SOLVER_NAME   = {"NONE", "GRIDDER", "LBM", "LPT"};
static constexpr std::array<std::string_view, 4> SOLVER_NAMELC = {"none", "gridder", "lbm", "lpt"};

static const std::vector<std::vector<GDouble>> DEFAULT_BOUNDINGBOX = {
    {0.0, 1.0}, {0.0, 1.0, 0.0, 1.0}, {0.0, 1.0, 0.0, 1.0, 0.0, 1.0}, {0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0}};

// just some spaces to arrange output
static constexpr std::string_view SP1{"  "};
static constexpr std::string_view SP2{"    "};
static constexpr std::string_view SP3{"      "};
static constexpr std::string_view SP4{"        "};
static constexpr std::string_view SP5{"          "};
static constexpr std::string_view SP6{"            "};
static constexpr std::string_view SP7{"              "};

namespace coordinate {
static constexpr std::array<char, 4> name = {'x', 'y', 'z', 'w'};
}


// Note if you adjust the GeomTypes also adjust GeomTypeString and resolveGeomType()!!!!
enum class GeomType { sphere, cube, box, stl, unknown, NumTypes };
static constexpr std::array<std::string_view, static_cast<GInt>(GeomType::NumTypes)> GeomTypeString = {"sphere", "cube", "box", "stl",
                                                                                                       "unknown"};

static inline auto resolveGeomType(const GString& type) -> GeomType {
  GInt index = std::distance(GeomTypeString.begin(), std::find(GeomTypeString.begin(), GeomTypeString.end(), type));
  if(index == static_cast<GInt>(GeomType::sphere)) {
    return GeomType::sphere;
  }
  if(index == static_cast<GInt>(GeomType::cube)) {
    return GeomType::cube;
  }
  if(index == static_cast<GInt>(GeomType::box)) {
    return GeomType::box;
  }
  if(index == static_cast<GInt>(GeomType::stl)) {
    return GeomType::stl;
  }
  return GeomType::unknown;
}

enum class DirId { mX, pX, mY, pY, mZ, pZ };
static constexpr std::array<std::string_view, 6> DirIdString = {"-x", "+x", "-y", "+y", "-z", "+z"};

static inline auto dirIdString2Id(const GString& dirId) -> GInt {
  if(dirId == "-x") {
    return 0;
  }
  if(dirId == "+x") {
    return 1;
  }
  if(dirId == "-y") {
    return 2;
  }
  if(dirId == "+y") {
    return 3;
  }
  if(dirId == "-z") {
    return 4;
  }
  if(dirId == "+z") {
    return 5;
  }
  std::cerr << "ERROR: Invalid direction " << dirId << std::endl;
  std::exit(-1);
}

#endif
