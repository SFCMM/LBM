#ifndef GRIDGENERATOR_APP_INTERFACE_H
#define GRIDGENERATOR_APP_INTERFACE_H

#include "interface/grid_interface.h"

class SolverInterface {
 public:
  virtual ~SolverInterface() = default;

  virtual void init(int argc, GChar** argv, GString config_file) = 0;
  virtual void initBenchmark(int argc, GChar** argv)             = 0;

  virtual auto run() -> GInt                           = 0;
  virtual auto grid() const -> const GridInterface&    = 0;
  virtual void transferGrid(const GridInterface& grid) = 0;
};

#endif // GRIDGENERATOR_APP_INTERFACE_H
