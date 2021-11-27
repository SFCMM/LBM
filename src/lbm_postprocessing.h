#ifndef LBM_LBM_POSTPROCESSING_H
#define LBM_LBM_POSTPROCESSING_H
#include "postprocessing_func.h"

template <GInt NDIM>
class LBMPostprocessFunctionExecutor : public PostprocessFunctionExecutorInterface<NDIM>,
                                       private PostprocessCartesianFunctionExecutor<NDIM> {
  LBMPostprocessFunctionExecutor(const pp::FuncType type)
    : PostprocessCartesianFunctionExecutor<NDIM>(type){

    };
};

#endif // LBM_LBM_POSTPROCESSING_H
