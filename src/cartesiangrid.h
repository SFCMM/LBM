#ifndef GRIDGENERATOR_CARTESIANGRID_H
#define GRIDGENERATOR_CARTESIANGRID_H

#include <gcem.hpp>

#include <sfcmm_common.h>
#include "base_cartesiangrid.h"
#include "cartesiangrid_generation.h"
#include "celltree.h"
#include "geometry.h"
#include "globaltimers.h"
#include "gridInterface.h"
#include "IO.h"

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class CartesianGrid : public BaseCartesianGrid<DEBUG_LEVEL, NDIM> {
 public:
  CartesianGrid()                     = default;
  ~CartesianGrid() override           = default;
  CartesianGrid(const CartesianGrid&) = delete;
  CartesianGrid(CartesianGrid&&)      = delete;
  auto operator=(const CartesianGrid&) -> CartesianGrid& = delete;
  auto operator=(CartesianGrid&&) -> CartesianGrid& = delete;

  void setCapacity(const GInt capacity) override {
    if(!m_tree.empty()) {
      TERMM(-1, "Invalid operation tree already allocated.");
    }
    m_tree.reset(capacity);
  }

  /// Load the generated grid in-memory and set additional properties
  /// \param grid Generated grid.
  void loadGridInplace(const CartesianGridGen<DEBUG_LEVEL, NDIM>& grid) {
    // grid.balance(); //todo: implement
  }

 private:
  cartesian::Tree<DEBUG_LEVEL, NDIM> m_tree;
};

#endif // GRIDGENERATOR_CARTESIANGRID_H
