// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_INTEGRATION_H
#define SFCMM_INTEGRATION_H

#include "common/sfcmm_types.h"

namespace integration {

struct range {
  GDouble a;
  GDouble b;
};

inline auto midpoint(const range range, const std::function<GDouble(GDouble)>& f) -> GDouble {
  assert(range.a < range.b);
  const GDouble h        = range.b - range.a;
  const GDouble midpoint = range.a + 0.5 * h;
  return h * f(midpoint);
}

inline auto multi_midpoint(const range range, const std::function<GDouble(GDouble)>& f, const GInt num_splits = 2) -> GDouble {
  assert(range.a < range.b);
  const GDouble h      = (range.b - range.a) / static_cast<GDouble>(num_splits);
  GDouble       result = 0;
  for(GInt i = 0; i < num_splits; ++i) {
    result += h * f(range.a + 0.5 * h * static_cast<GDouble>(2 * i + 1));
  }
  return result;
}

inline auto trapezoidal(const range range, const std::function<GDouble(GDouble)>& f) -> GDouble {
  assert(range.a < range.b);
  const GDouble h  = range.b - range.a;
  const GDouble x1 = range.a;
  const GDouble x2 = range.b;
  return 0.5 * h * (f(x1) + f(x2));
}

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

inline auto simpson(const range range, const std::function<GDouble(GDouble)>& f) -> GDouble {
  assert(range.a < range.b);
  const GDouble h  = 0.5 * (range.b - range.a);
  const GDouble x1 = range.a;
  const GDouble x2 = range.a + h;
  const GDouble x3 = range.b;
  return 1.0 / 3.0 * h * (f(x1) + 4 * f(x2) + f(x3));
}

inline auto multi_simpson(const range range, const std::function<GDouble(GDouble)>& f, const GInt num_splits = 2) -> GDouble {
  assert(range.a < range.b);
  const GDouble h      = 0.5 * (range.b - range.a) / static_cast<GDouble>(num_splits);
  const GDouble x1     = range.a;
  const GDouble x2     = range.a + h;
  const GDouble x3     = range.b;
  GDouble       result = 0;
  for(GInt i = 0; i < num_splits; ++i) {
    result += 1.0 / 3.0 * h * (f(range.a + 2 * i * h) + 4 * f(range.a + (2 * i + 1) * h) + f(range.a + (2 * i + 2) * h));
  }
  return result;
}

} // namespace integration

#endif // SFCMM_INTEGRATION_H
