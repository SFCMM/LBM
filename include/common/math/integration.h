// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_INTEGRATION_H
#define SFCMM_INTEGRATION_H

#include "common/sfcmm_types.h"

namespace integration {

struct range {
  GDouble a;
  GDouble b;
};

/// Integrate given function by using the midpoint rule (1 Function eval) Accuracy O(N)
/// \param range Integration range
/// \param f function to integrate
/// \return Integrated value
inline auto midpoint(const range range, const std::function<GDouble(GDouble)>& f) -> GDouble {
  assert(range.a < range.b);
  const GDouble h        = range.b - range.a;
  const GDouble midpoint = range.a + 0.5 * h;
  return h * f(midpoint);
}

/// Integrate by splitting range into multiple equidistant intervals. (num_splits function evals) Accuracy O(N)
/// \param range Integration range
/// \param f function to integrate
/// \param num_splits Split the given integration range into num_split intervals
/// \return Integrated value
inline auto multi_midpoint(const range range, const std::function<GDouble(GDouble)>& f, const GInt num_splits = 2) -> GDouble {
  assert(range.a < range.b);
  const GDouble h      = (range.b - range.a) / static_cast<GDouble>(num_splits);
  GDouble       result = 0;
  for(GInt i = 0; i < num_splits; ++i) {
    result += h * f(range.a + 0.5 * h * static_cast<GDouble>(2 * i + 1));
  }
  return result;
}

/// Integrate given function by using the trapezoidal rule (2 Function evals) Accuracy O(N)
/// \param range Integration range
/// \param f function to integrate
/// \return Integrated value
inline auto trapezoidal(const range range, const std::function<GDouble(GDouble)>& f) -> GDouble {
  assert(range.a < range.b);
  const GDouble h  = range.b - range.a;
  const GDouble x1 = range.a;
  const GDouble x2 = range.b;
  return 0.5 * h * (f(x1) + f(x2));
}

/// Integrate by splitting range into multiple equidistant intervals. (num_splits + 1 function evals) Accuracy O(N)
/// \param range Integration range
/// \param f function to integrate
/// \param num_splits Split the given integration range into num_split intervals
/// \return Integrated value
inline auto multi_trapezoidal(const range range, const std::function<GDouble(GDouble)>& f, const GInt num_splits = 2) -> GDouble {
  assert(range.a < range.b);
  const GDouble h      = (range.b - range.a) / static_cast<GDouble>(num_splits);
  GDouble       fprev  = f(range.a);
  GDouble       result = 0;
  for(GInt i = 0; i < num_splits; ++i) {
    result += 0.5 * h * fprev;
    fprev = f(range.a + h * (i + 1));
    result += 0.5 * h * fprev;
  }
  return result;
}

/// Integrate given function by using the simpson rule (3 Function evals) Accuracy O(N^3)
/// \param range Integration range
/// \param f function to integrate
/// \return Integrated value
inline auto simpson(const range range, const std::function<GDouble(GDouble)>& f) -> GDouble {
  assert(range.a < range.b);
  const GDouble h  = 0.5 * (range.b - range.a);
  const GDouble x1 = range.a;
  const GDouble x2 = range.a + h;
  const GDouble x3 = range.b;
  return 1.0 / 3.0 * h * (f(x1) + 4 * f(x2) + f(x3));
}

/// Integrate by splitting range into multiple equidistant intervals. (2*num_splits+1 function evals) Accuracy O(N^3)
/// \param range Integration range
/// \param f function to integrate
/// \param num_splits Split the given integration range into num_split intervals
/// \return Integrated value
inline auto multi_simpson(const range range, const std::function<GDouble(GDouble)>& f, const GInt num_splits = 2) -> GDouble {
  assert(range.a < range.b);
  const GDouble h      = 0.5 * (range.b - range.a) / static_cast<GDouble>(num_splits);
  GDouble       result = 0;
  GDouble       f1     = 0;
  GDouble       f2     = 0;
  GDouble       f3     = f(range.a);
  for(GInt i = 0; i < num_splits; ++i) {
    f1 = f3;
    f2 = f(range.a + (2 * i + 1) * h);
    f3 = f(range.a + (2 * i + 2) * h);
    result += 1.0 / 3.0 * h * (f1 + 4 * f2 + f3);
  }
  return result;
}

} // namespace integration

#endif // SFCMM_INTEGRATION_H
