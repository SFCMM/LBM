// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_CARTESIAN_H
#define SFCMM_CARTESIAN_H

#include <gcem.hpp>
#include <cassert>

namespace cartesian {

template <GInt NDIM>
static constexpr inline auto maxNoNghbrsDiag() -> GInt;

/// The direction opposite to the given dir. [Main directions only]
/// \param dir Direction for which to obtain opposite direction.
/// \return The opposite direction to dir.
static constexpr inline auto oppositeDir(const GInt dir) -> GInt { return dir + 1 - 2 * (dir % 2); }

/// The direction opposite to the given dir. [Including diagonals]
/// \param dir Direction for which to obtain opposite direction.
/// \return The opposite direction to dir.
template <GInt NDIM>
static constexpr inline auto oppositeDir(const GInt dir) -> GInt {
  assert(dir < 2 * maxNoNghbrsDiag<NDIM>() && "Invalid direction");
  if(dir < 2 * NDIM) {
    return dir + 1 - 2 * (dir % 2);
  }
  if constexpr(NDIM == 2) {
    // diagonal direction
    return dir > 5 ? dir - 2 : dir + 2;
  }
  std::cerr << "Invalid dir in oppositeDir() " << std::endl;
  std::exit(-1);
}

/// Return maximum number of children per cell
/// \tparam NDIM Number of dimension of the Cartesian grid
/// \return The maximum number of children for a NDIM Cartesian mesh
template <GInt NDIM>
static constexpr inline auto maxNoChildren() -> GInt {
  return gcem::pow(2, NDIM);
}

/// Return maximum number of Neighbors per cell
/// \tparam NDIM Number of dimensions of the Cartesian grid
/// \return The maximum number of same level neighbors at the given NDIM.
template <GInt NDIM>
static constexpr inline auto maxNoNghbrs() -> GInt {
  return 2 * NDIM;
}

/// Return maximum number of Neighbors including diagonals (2D/3D) and tridiagonals (3D) per cell
/// \tparam NDIM Number of dimensions of the Cartesian grid
/// \return The maximum number of same level neighbors at the given NDIM including diagonals.
template <GInt NDIM>
static constexpr inline auto maxNoNghbrsDiag() -> GInt {
  if(NDIM == 1) {
    return 2;
  }
  if(NDIM == 2) {
    return 8;
  }
  if(NDIM == 3) {
    return 26;
  }
  return 50;
}

// todo: replace with constant expression function
/// Given the childId gives the "direction" of this child relative to the center
/// of a cell.
static constexpr std::array<std::array<GDouble, sfcmm::MAX_DIM>, cartesian::maxNoChildren<sfcmm::MAX_DIM>()> childDir = {{
    //-> 2D
    // -x,-y, -z
    {{-1, -1, -1, -1}}, // 0
    //+x, -y, -z
    {{1, -1, -1, -1}}, // 1
    //-x, +y, -z
    {{-1, 1, -1, -1}}, // 2
    //+x, -y, -z
    {{1, 1, -1, -1}}, // 3
    //<- 2D

    //-> 3D (+z)
    {{-1, -1, 1, -1}}, // 4
    {{1, -1, 1, -1}},  // 5
    {{-1, 1, 1, -1}},  // 6
    {{1, 1, 1, -1}},   // 7
    //<- 3D (+z)

    //-> 4D
    {{-1, -1, -1, 1}}, // 8
    {{1, -1, -1, 1}},  // 9
    {{-1, 1, -1, 1}},  // 10
    {{1, 1, -1, 1}},   // 11
    {{-1, -1, 1, 1}},  // 12
    {{1, -1, 1, 1}},   // 13
    {{-1, 1, 1, 1}},   // 14
    {{1, 1, 1, 1}}     // 15
                       //<- 4D
}};

// todo: replace with constant expression function
/// Given the childId gives the neighboring childIds(and existence ==-1 ->
/// doesnot exist)
static constexpr std::array<std::array<GDouble, cartesian::maxNoNghbrs<sfcmm::MAX_DIM>()>, cartesian::maxNoChildren<sfcmm::MAX_DIM>()>
    nghbrInside = {{
        //-x +x -y +y -z +z -zz +zz
        {{-1, 1, -1, 2, -1, 4, -1, 8}},  // 0
        {{0, -1, -1, 3, -1, 5, -1, 9}},  // 1
        {{-1, 3, 0, -1, -1, 6, -1, 10}}, // 2
        {{2, -1, 1, -1, -1, 7, -1, 11}}, // 3
        {{-1, 5, -1, 6, 0, -1, -1, 12}}, // 4
        {{4, -1, -1, 7, 1, -1, -1, 13}}, // 5
        {{-1, 7, 4, -1, 2, -1, -1, 14}}, // 6
        {{6, -1, 5, -1, 3, -1, -1, 15}}, // 7

        // upper table +8 and last dir = -1
        {{-1, 9, -1, 10, -1, 12, 0, -1}},  // 8
        {{8, -1, -1, 11, -1, 13, 1, -1}},  // 9
        {{-1, 11, 8, -1, -1, 14, 2, -1}},  // 10
        {{10, -1, 9, -1, -1, 15, 3, -1}},  // 11
        {{-1, 13, -1, 14, 8, -1, 4, -1}},  // 12
        {{12, -1, -1, 15, 9, -1, 5, -1}},  // 13
        {{-1, 15, 12, -1, 10, -1, 6, -1}}, // 14
        {{14, -1, 13, -1, 11, -1, 7, -1}}  // 15
    }};

// todo: replace with constant expression function
/// Given the childId obtain the possible neighbors in a neighboring cell that
/// doesnot have the same parent
static constexpr std::array<std::array<GDouble, cartesian::maxNoNghbrs<sfcmm::MAX_DIM>()>, cartesian::maxNoChildren<sfcmm::MAX_DIM>()>
    nghbrParentChildId = {{
        //-x +x -y +y -z +z -zz +zz
        {{1, -1, 2, -1, 4, -1, 8, -1}},  // 0
        {{-1, 0, 3, -1, 5, -1, 9, -1}},  // 1
        {{3, -1, -1, 0, 6, -1, 10, -1}}, // 2
        {{-1, 2, -1, 1, 7, -1, 11, -1}}, // 3
        {{5, -1, 6, -1, -1, 0, 12, -1}}, // 4
        {{-1, 4, 7, -1, -1, 1, 13, -1}}, // 5
        {{7, -1, -1, 4, -1, 2, 14, -1}}, // 6
        {{-1, 6, -1, 5, -1, 3, 15, -1}}, // 7

        // upper table +8 and last dir = -1
        {{9, -1, 10, -1, 12, -1, -1, 0}},  // 8
        {{-1, 8, 11, -1, 13, -1, -1, 1}},  // 9
        {{11, -1, -1, 8, 14, -1, -1, 2}},  // 10
        {{-1, 10, -1, 9, 15, -1, -1, 3}},  // 11
        {{13, -1, 14, -1, -1, 8, -1, 4}},  // 12
        {{-1, 12, 15, -1, -1, 9, -1, 5}},  // 13
        {{15, -1, -1, 12, -1, 10, -1, 6}}, // 14
        {{-1, 14, -1, 13, -1, 11, -1, 7}}  // 15
    }};

/// Return the direction diagonal between 2 (2D) directions
/// \tparam NDIM Dimensionality
/// \param dir1 Direction 1
/// \param dir2 Direction 2
/// \return Direction diagonal between these
template <GInt NDIM>
static constexpr auto inbetweenDiagDirs(GInt /*dir1*/, GInt /*dir2*/) -> GInt {
  std::cerr << "Invalid dist in inbetweenDiagDirs() " << std::endl;
  std::exit(-1);
  return -1;
}

template <>
constexpr auto inbetweenDiagDirs<2>(GInt dir1, GInt dir2) -> GInt {
  if(dir2 < dir1) {
    GInt tmp = dir2;
    dir2     = dir1;
    dir1     = tmp;
  }
  // after this dist 1  < dir2
  switch(dir1) {
    case 0:
      if(dir2 == 2) {
        return 6;
      }
      if(dir2 == 3) {
        return 7;
      }
    case 1:
      if(dir2 == 2) {
        return 5;
      }
      if(dir2 == 3) {
        return 4;
      }
      break;
    default:
      std::cerr << "Invalid dist in inbetweenDiagDirs<2>()" << std::endl;
      std::exit(-1);
  }
  return -1;
}

/// Return the direction unit vector.
/// \tparam NDIM Dimensionality
/// \param dir Direction
/// \return Unit vector of the corresponding direction
template <GInt NDIM>
inline auto dirVec(GInt /*dir*/) -> VectorD<NDIM> {
  std::cerr << "Invalid dir in dirVec()" << std::endl;
  std::exit(-1);
  VectorD<NDIM> tmp;
  tmp.fill(NAN);
  return tmp;
}

template <>
inline auto dirVec<1>(GInt dir) -> VectorD<1> {
  VectorD<1> tmp;
  switch(dir) {
    case 0:
      tmp.fill(-1);
      return tmp;
    case 1:
      tmp.fill(1);
      return tmp;
    default:
      std::cerr << "Invalid dir in dirVec<1>()" << std::endl;
      std::exit(-1);
  }
  tmp.fill(NAN);
  return tmp;
}

template <>
inline auto dirVec<2>(GInt dir) -> VectorD<2> {
  switch(dir) {
    case 0:
      return {-1, 0};
    case 1:
      return {1, 0};
    case 2:
      return {0, -1};
    case 3:
      return {0, 1};
    case 4:
      return {1, 1};
    case 5:
      return {1, -1};
    case 6:
      return {-1, -1};
    case 7:
      return {-1, 1};
    default:
      std::cerr << "Invalid dir in dirVec<2>()" << std::endl;
      std::exit(-1);
  }
  return {NAN, NAN};
}

} // namespace cartesian
#endif // SFCMM_CARTESIAN_H
