#ifndef LBM_LBM_SOLVER_H
#define LBM_LBM_SOLVER_H
#include "interface/app_interface.h"

template <Debug_Level DEBUG_LEVEL>
class LBMSolver : public AppInterface {
 public:
  void               init(int argc, GChar** argv, GString config_file) override{};
  void               initBenchmark(int argc, GChar** argv) override{};
  auto               run() -> GInt override{};
  [[nodiscard]] auto grid() const -> const GridInterface& override { return *m_grid; };
  void               transferGrid(const GridInterface& grid) const override {
    switch(m_dim) {
      case 1:
        static_cast<CartesianGrid<DEBUG_LEVEL, 1>*>(m_grid.get())
            ->loadGridInplace(*static_cast<const CartesianGridGen<DEBUG_LEVEL, 1>*>(static_cast<const void*>(&grid)));
        break;
      case 2:
        static_cast<CartesianGrid<DEBUG_LEVEL, 2>*>(m_grid.get())
            ->loadGridInplace(*static_cast<const CartesianGridGen<DEBUG_LEVEL, 2>*>(static_cast<const void*>(&grid)));
        break;
      case 3:
        static_cast<CartesianGrid<DEBUG_LEVEL, 3>*>(m_grid.get())
            ->loadGridInplace(*static_cast<const CartesianGridGen<DEBUG_LEVEL, 3>*>(static_cast<const void*>(&grid)));
        break;
      default:
        TERMM(-1, "Only dimensions 1,2 and 3 are supported.");
    }
  };

 private:
  std::unique_ptr<GridInterface> m_grid;

  GInt m_dim = 0;
};

#endif // LBM_LBM_SOLVER_H
