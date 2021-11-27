#ifndef LBM_POSTPROCESSING_CARTESIAN_H
#define LBM_POSTPROCESSING_CARTESIAN_H
#include "postprocessing_func.h"

template<GInt NDIM>
class PostprocessCartesianFunctionExecutor{
 public:
  PostprocessCartesianFunctionExecutor(const pp::FuncType type): m_type(type){

  };

  void setup(const pp::FuncType type){
    m_type = type;

  }

  void execute(){
    if(m_type == pp::FuncType::LINE){
      line();
    }
  }

 private:
  PPFuncType m_type;

};

#endif // LBM_POSTPROCESSING_CARTESIAN_H
