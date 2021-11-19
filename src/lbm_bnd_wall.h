#ifndef LBM_BND_WALL_H
#define LBM_BND_WALL_H
#include "lbm_bnd_interface.h"
#include "lbm_constants.h"

template <LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBndCell_wallBB : public LBMBndCell<LBTYPE> {
  // class LBMBndCell_wallBB {
 public:
  LBMBndCell_wallBB() = delete;
  //  LBMBndCell_wallBB() = default;
  LBMBndCell_wallBB(const GInt cellId, const VectorD<dim(LBTYPE)>& _normal) : LBMBndCell<LBTYPE>(cellId, _normal) {}
  //  LBMBndCell_wallBB(const GInt cellId) {}
  virtual ~LBMBndCell_wallBB() = default;

  void init() override {
    LBMBndCell<LBTYPE>::init();

    // precalculate weight and bndIndex
    for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
      // determine if the dist points into the outside normal direction of bndry
      if(inDirection<dim(LBTYPE)>(normal(), LBMethod<LBTYPE>::m_dirs[dist])) {
        m_bndIndex[m_noSetDists] = dist;
        ++m_noSetDists;
      }
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

      if constexpr(TANGENTIALVELO) {
        fold(mapped(), oppositeDist) += m_tangentialVelo[oppositeDist];
      }
    }
  }

  void setTangentialVelocity(const GDouble tangentialVelo) {
    if constexpr(TANGENTIALVELO) {
      static constexpr GDouble parallelRequirement = 10 * GDoubleEps;

      m_tangentialVelo.fill(0);
      const auto& _normal = normal();
      for(GInt id = 0; id < m_noSetDists; ++id) {
        const GInt           dist = LBMethod<LBTYPE>::oppositeDist(m_bndIndex[id]);
        VectorD<dim(LBTYPE)> dir;
        for(GInt n = 0; n < dim(LBTYPE); ++n) {
          dir[n] = LBMethod<LBTYPE>::m_dirs[dist][n];
        }

        VectorD<dim(LBTYPE)> tangentialVec;
        if(dim(LBTYPE) == 2) {
          tangentialVec[0] = _normal[1];
          tangentialVec[1] = _normal[0];
        } else {
          TERMM(-1, "Not implemented");
        }


        const GDouble normalDotDir     = _normal.dot(dir);
        const GDouble tangentialDotDir = tangentialVec.dot(dir);
        const GDouble angle            = gcem::acos(normalDotDir / dir.norm());
        const GBool   parallel         = std::abs(angle - PI) < parallelRequirement;

        // vector is parallel to the normal vector -> velocity = 0
        if(!parallel) {
          m_tangentialVelo[dist] = tangentialDotDir * 1.0 / 6.0 * tangentialVelo;
        }
      }
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

  std::array<GDouble, noDists(LBTYPE) * static_cast<GInt>(TANGENTIALVELO)> m_tangentialVelo;
  std::array<GInt, noDists(LBTYPE)>                                        m_bndIndex;
  GInt                                                                     m_noSetDists = 0;
};

template <LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBnd_wallBB : public LBMBndInterface {
 public:
  LBMBnd_wallBB(const Surface<dim(LBTYPE)>& surf, const json& properties) {
    GInt surfId = 0;
    for(const GInt cellId : surf.getCellList()) {
      m_bndCells.emplace_back(cellId, surf.normal(surfId));
      ++surfId;
    }
    if constexpr(TANGENTIALVELO) {
      m_tangentialVelo = config::required_config_value<GDouble>(properties, "tangentialVelocity");
      logger << "Setting tangentialVelocity " << m_tangentialVelo << std::endl;
    }
    LBMBnd_wallBB<LBTYPE, TANGENTIALVELO>::init();
  }
  ~LBMBnd_wallBB() override = default;

  // deleted constructors not needed
  LBMBnd_wallBB(const LBMBnd_wallBB&) = delete;
  LBMBnd_wallBB(LBMBnd_wallBB&&)      = delete;
  auto operator=(const LBMBnd_wallBB&) -> LBMBnd_wallBB& = delete;
  auto operator=(LBMBnd_wallBB&&) -> LBMBnd_wallBB& = delete;

  void init() override {
    for(auto& bndCell : m_bndCells) {
      bndCell.init();
      bndCell.setTangentialVelocity(m_tangentialVelo);
    }
  }

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold, const std::function<GDouble&(GInt, GInt)>& vars) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& vars) override {
    for(auto& bndCell : m_bndCells) {
      bndCell.apply(f, fold);
    }
  }

 private:
  GDouble                                                m_tangentialVelo = 0.1;
  std::vector<LBMBndCell_wallBB<LBTYPE, TANGENTIALVELO>> m_bndCells;
};
#endif // LBM_BND_WALL_H
