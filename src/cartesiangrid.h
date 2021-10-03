#ifndef GRIDGENERATOR_CARTESIANGRID_H
#define GRIDGENERATOR_CARTESIANGRID_H

#include <gcem.hpp>

#include <sfcmm_common.h>
#include "base_cartesiangrid.h"
#include "celltree.h"
#include "common/IO.h"
#include "geometry.h"
#include "globaltimers.h"
#include "gridgenerator/cartesiangrid_generation.h"
#include "interface/grid_interface.h"

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
    m_tree.setCapacity(capacity);
  }

  void reset() override { m_tree.reset(); }

  void save(const GString& fileName, const json& gridOutConfig) override { TERMM(-1, "Not implemented!"); }

  /// Load the generated grid in-memory and set additional properties
  /// \param grid Generated grid.
  void loadGridInplace(const CartesianGridGen<DEBUG_LEVEL, NDIM>& grid) {
    // grid.balance(); //todo: implement
    m_tree.setCapacity(grid.capacity()); // todo: change for adaptation

#ifdef _OPENMP
#pragma omp parallel default(none) shared(grid)
    {
#endif
      const GInt noCells = grid.size();
#ifdef _OPENMP
#pragma omp for
#endif
      for(GInt cellId = 0; cellId < noCells; ++cellId) {
        m_tree.globalId(cellId) = grid.globalId(cellId);
        m_tree.parent(cellId)   = grid.parent(cellId);
        for(GInt childId = 0; childId < cartesian::maxNoChildren<NDIM>(); ++childId) {
          m_tree.child(cellId, childId) = grid.child(cellId, childId);
        }
        for(GInt nghbrId = 0; nghbrId < cartesian::maxNoNghbrs<NDIM>(); ++nghbrId) {
          m_tree.neighbor(cellId, nghbrId) = grid.neighbor(cellId, nghbrId);
        }
        m_tree.level(cellId) = grid.level(cellId);
        for(GInt dir = 0; dir < NDIM; ++dir) {
          m_tree.center(cellId, dir) = grid.center(cellId, dir);
        }
      }
#ifdef _OPENMP
    }
#endif
    setProperties();
    if(m_loadBalancing) {
      setWorkload();
      calculateOffspringsAndWeights();
    }
  }

 private:
  void setProperties(){

  };
  void setWorkload() { TERMM(-1, "Not implemented!"); };
  void calculateOffspringsAndWeights() { TERMM(-1, "Not implemented!"); };

  cartesian::Tree<DEBUG_LEVEL, NDIM> m_tree{};

  GBool m_loadBalancing = false;
};

#endif // GRIDGENERATOR_CARTESIANGRID_H
