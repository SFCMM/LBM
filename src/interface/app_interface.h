#ifndef GRIDGENERATOR_APP_INTERFACE_H
#define GRIDGENERATOR_APP_INTERFACE_H

class AppInterface {
 public:
  virtual void init(int argc, GChar** argv, GString config_file) = 0;
  virtual void initBenchmark(int argc, GChar** argv)             = 0;
  virtual auto run() -> GInt                                     = 0;
  virtual auto grid() const -> const GridInterface&              = 0;
  virtual void transferGrid(const GridInterface& grid) const     = 0;
};

#endif // GRIDGENERATOR_APP_INTERFACE_H
