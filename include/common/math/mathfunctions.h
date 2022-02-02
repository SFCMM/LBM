// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_MATHFUNCTIONS_H
#define SFCMM_MATHFUNCTIONS_H
#include "common/sfcmm_types.h"
#include <cmath>
#include <gcem.hpp>

namespace detail_ {
template <class T, class U> struct APPROX_ERROR {};
} // namespace detail_

/// Approximately equality between two floating numbers.
/// \tparam T
/// \tparam U
/// \return
template <class T, class U>
constexpr inline auto approx(const T& /*unused*/, const U& /*unused*/, const T& /*unused*/) -> bool {
  using error = typename detail_::APPROX_ERROR<T, U>::ERROR_BOTH_TYPES_MUST_BE_GDouble_or_GFloat;
  error();
  return true;
}

/// Approximately equality between two doubles.
/// \param a First double.
/// \param b Second double.
/// \param eps Precision used for comparison.
/// \return If a is equal to b within eps.
template <>
constexpr inline auto approx<double, double>(const double& a, const double& b, const double& eps) -> bool {
  return gcem::abs(a - b) < eps;
}

/// Approximately equality between two floats.
/// \param a First float.
/// \param b Second float.
/// \param eps Precision used for comparison.
/// \return If a is equal to b within eps.
template <>
constexpr inline auto approx<float, float>(const float& a, const float& b, const float& eps) -> bool {
  return gcem::abs(a - b) < eps;
}

/// Approximately equality between two floating numbers.
/// \tparam T
/// \tparam U
/// \return
template <class T, class U>
constexpr inline auto approx(const T& /*unused*/, const U& /*unused*/) -> bool {
  using error = typename detail_::APPROX_ERROR<T, U>::ERROR_BOTH_TYPES_MUST_BE_GDouble_or_GFloat;
  error();
  return true;
}

/// Approximately equality between two doubles.
/// \param a First double.
/// \param b Second double.
/// \return If a is equal to b within eps.
template <>
constexpr inline auto approx<double, double>(const double& a, const double& b) -> bool {
  return gcem::abs(a - b) < std::numeric_limits<double>::epsilon();
}

/// Approximately equality between two floats.
/// \param a First float.
/// \param b Second float.
/// \return If a is equal to b within eps.
template <>
constexpr inline auto approx<float, float>(const float& a, const float& b) -> bool {
  return gcem::abs(a - b) < std::numeric_limits<double>::epsilon();
}

/// Check if a value is even.
/// \tparam T Integer type
/// \param num Number to be checked to be even
/// \return True if it is an even number
template <class T>
[[nodiscard]] static constexpr inline auto isEven(T num) -> bool {
  return num % 2 == 0;
}

/// Convert degrees to radians.
/// \param degrees Angle in degrees
/// \return Angle in Radians
[[nodiscard]] static constexpr inline auto toRadians(const GDouble degrees) { return degrees / 180.0 * M_PI; }


#endif // SFCMM_MATHFUNCTIONS_H
