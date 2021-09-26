#ifndef LBM_LBM_SOLVER_H
#define LBM_LBM_SOLVER_H
#include "app_interface.h"

template <Debug_Level DEBUG_LEVEL>
class LBMSolver : public AppInterface {
 public:
  void init(int argc, GChar** argv, GString config_file) {};
  void initBenchmark(int argc, GChar** argv)             {};
  auto run() -> GInt                                     {};
};

#endif // LBM_LBM_SOLVER_H
