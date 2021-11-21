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

    // set which dists are to be set by this boundary condition
    for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
      // determine if the dist points into the outside normal direction of bndry
      if(inDirection<dim(LBTYPE)>(normal(), LBMethod<LBTYPE>::m_dirs[dist])) {
        m_bndIndex[m_noSetDists] = dist;
        ++m_noSetDists;
      }
    }

    //todo: this needs to be moved since we cannot guarantee double matching cells are handled correctly
    for(GInt id = 0; id < m_noSetDists; ++id) {
      const GInt dist = m_bndIndex[id];
      auto centerA        = surfConnected->center(mapped());
      for(GInt dir = 0; dir < dim(LBTYPE); ++dir){
        centerA[dir] += LBMethod<LBTYPE>::m_dirs[dist][dir] * 0.5 * surfConnected->cellLength(mapped());
      }


      GInt minDistCellId = -1;
      GDouble minDist = std::numeric_limits<GDouble>::max();
      // connect cell of the boundary surface and connected surface
      for(const GInt cellIdB : surfConnected->getCellList()) {
        const auto& centerB        = surfConnected->center(cellIdB);
        const GDouble diff = (centerA - centerB).norm();

        if(minDist > diff){
          minDist = diff;
          minDistCellId = cellIdB;
        }
      }

      m_linkedCell[id] = minDistCellId;
      ASSERT(minDistCellId >= 0, "Invalid Cellid");
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) {
    // iterate over the distributions that need to be set
    for(GInt id = 0; id < m_noSetDists; ++id) {
      // dists that point into the outside direction
      const GInt dist = m_bndIndex[id];
      // dist in reflected/opposite direction i.e. to the inside of the wall
      GInt oppositeDist = LBMethod<LBTYPE>::oppositeDist(dist);

      // standard bounceback i.e. distribution that hits the wall is reflected to the opposite distribution direction
      fold(mapped(), oppositeDist) = f(mapped(), dist);
    }
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
                const std::function<GDouble&(GInt, GInt)>& vars) override {}
  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& vars) override {}

 private:
  std::vector<LBMBndCell_periodic<LBTYPE>> m_bndCells;
};
#endif // LBM_LBM_BND_PERIODIC_H
