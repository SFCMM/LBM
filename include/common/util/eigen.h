// SPDX-License-Identifier: BSD-3-Clause

#ifndef COMMON_EIGEN_H
#define COMMON_EIGEN_H
#include "common/sfcmm_types.h"
namespace eigenutil {
/// Unpack a std::vector<double> to an eigen compatible type
/// \tparam NDIM Dimensions of the eigen type
/// \param vec Vector to be unpacked
/// \return Eigen compatible vector type
template <GInt NDIM>
static constexpr auto unpack(const std::vector<GDouble>& vec) -> VectorD<NDIM> {
  VectorD<NDIM> tmpVec;
  for(GUint i = 0; i < vec.size(); ++i) {
    tmpVec[i] = vec[i];
  }
  return tmpVec;
}
} // namespace eigenutil

#endif // COMMON_EIGEN_H
