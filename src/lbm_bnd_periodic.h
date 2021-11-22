#ifndef LBM_LBM_BND_PERIODIC_H
#define LBM_LBM_BND_PERIODIC_H
#include "lbm_bnd_interface.h"
#include "lbm_constants.h"
#include "pv.h"

template <LBMethodType LBTYPE>
class LBMBndCell_periodic : public LBMBndCell<LBTYPE> {
  // class LBMBndCell_wallBB {
 public:
  LBMBndCell_periodic() = delete;
  LBMBndCell_periodic(const GInt cellId, const VectorD<dim(LBTYPE)>& _normal) : LBMBndCell<LBTYPE>(cellId, _normal) {}
  virtual ~LBMBndCell_periodic() = default;

  // todo: move to some common place
  void init(const Surface<dim(LBTYPE)>* surfConnected)  {
    LBMBndCell<LBTYPE>::init();

    const GDouble cellLength = surfConnected->cellLength(mapped());
    const auto center = surfConnected->center(mapped());
    const auto& bb = surfConnected->grid()->boundingBox();

    // set which dists are to be set by this boundary condition
    for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
      // determine if the dist points into the outside normal direction of bndry
      if(inDirection<dim(LBTYPE)>(normal(), LBMethod<LBTYPE>::m_dirs[dist])) {
        GBool inside = true;
        for(GInt dir = 0; dir < dim(LBTYPE); ++dir){
          if(std::abs(normal()[dir]) > 0){
            continue;
          }

          const GDouble coord =center[dir] + LBMethod<LBTYPE>::m_dirs[dist][dir] * cellLength;
          //todo: move to function isInside...
          if( coord < bb.min(dir) || coord > bb.max(dir)){
            inside = false;
            break;
          }
        }

        if(inside) {
          m_bndIndex[m_noSetDists] = dist;
          ++m_noSetDists;
        }
      }
    }

    //todo: this needs to be moved since we cannot guarantee double matching cells are handled correctly
    for(GInt id = 0; id < m_noSetDists; ++id) {
      m_linkedCell[id] = -1;
      const GInt dist = m_bndIndex[id];
      auto centerA        = center;
      for(GInt dir = 0; dir < dim(LBTYPE); ++dir){
        if(std::abs(normal()[dir]) > 0){
          continue;
        }
        centerA[dir] += LBMethod<LBTYPE>::m_dirs[dist][dir] * cellLength;
      }

      GDouble minDist = std::numeric_limits<GDouble>::max();
      // connect cell of the boundary surface and connected surface
      for(const GInt cellIdB : surfConnected->getCellList()) {
        const auto&   centerB = surfConnected->center(cellIdB);
        for(GInt dir = 0; dir < dim(LBTYPE); ++dir) {
          const GDouble diff    = std::abs(centerA[dir] - centerB[dir]);
          if(diff <= GDoubleEps) {
            m_linkedCell[id] = cellIdB;
            break;
          }
        }
        if(m_linkedCell[id] >=0){
          break;
        }
      }

      ASSERT(m_linkedCell[id] >= 0, "Invalid Cellid");
    }

    for(GInt idA = 0; idA < m_noSetDists; ++idA) {
      for(GInt idB = idA+1; idB < m_noSetDists; ++idB) {
        if(m_linkedCell[idA] == m_linkedCell[idB]){
          TERMM(-1, "Invalid periodic bnd");
        }
      }
    }
  }

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) {
    //do the same as in the propagation step
    for(GInt id = 0; id < m_noSetDists; ++id) {
      const GInt dist = m_bndIndex[id];
      fold(m_linkedCell[id], dist) = f(mapped(), dist);
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& vars) {
  }

  // todo: fix me
  //  deleted constructors not needed
  //  LBMBndCell_wallBB(const LBMBndCell_wallBB&) = delete;
  //  LBMBndCell_wallBB(LBMBndCell_wallBB&&)      = delete;
  //  auto operator=(const LBMBndCell_wallBB&) -> LBMBndCell_wallBB& = delete;
  //  auto operator=(LBMBndCell_wallBB&&) -> LBMBndCell_wallBB& = delete;

 private:
  using LBMBndCell<LBTYPE>::mapped;
  using LBMBndCell<LBTYPE>::normal;

  std::array<GInt, noDists(LBTYPE)> m_bndIndex{};
  std::array<GInt, noDists(LBTYPE)> m_linkedCell{};
  GInt                              m_noSetDists = 0;
};

template <LBMethodType LBTYPE>
class LBMBnd_Periodic : public LBMBndInterface {
 public:
  LBMBnd_Periodic(const Surface<dim(LBTYPE)>* surf, const Surface<dim(LBTYPE)>* surfConnected, const json& properties) {
    GInt surfId = 0;
    for(const GInt cellId : surf->getCellList()) {
      m_bndCells.emplace_back(cellId, surf->normal(surfId));
      ++surfId;
    }
    LBMBnd_Periodic<LBTYPE>::init(surfConnected);
  }

  void init(const Surface<dim(LBTYPE)>* surfConnected) {
    for(auto& bndCell : m_bndCells) {
      bndCell.init(surfConnected);
    }
  }

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                const std::function<GDouble&(GInt, GInt)>& vars) override {
    //apply to all boundary cells
    for(auto& bndCell : m_bndCells) {
      bndCell.preApply(f, fold);
    }
  }
  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& vars) override {}

 private:
  std::vector<LBMBndCell_periodic<LBTYPE>> m_bndCells;
};
#endif // LBM_LBM_BND_PERIODIC_H
