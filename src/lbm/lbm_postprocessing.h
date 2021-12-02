#ifndef LBM_LBM_POSTPROCESSING_H
#define LBM_LBM_POSTPROCESSING_H
#include "postprocessing_cartesian.h"

template <GInt NDIM>
class LBMPostprocessFunctionLine : public PostprocessCartesianFunctionLine<NDIM> {
 public:
  LBMPostprocessFunctionLine(const json& conf, const CartesianGridData<NDIM>& data) : PostprocessCartesianFunctionLine<NDIM>(conf, data) {}
};

#endif // LBM_LBM_POSTPROCESSING_H
