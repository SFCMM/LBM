#ifndef LBM_LBM_BND_H
#define LBM_LBM_BND_H

#include <sfcmm_common.h>
#include "common/surface.h"
#include "configuration.h"


class LBMBndInterface {
 public:
  LBMBndInterface()          = default;
  virtual ~LBMBndInterface() = default;

  // deleted constructors not needed
  LBMBndInterface(const LBMBndInterface&) = delete;
  LBMBndInterface(LBMBndInterface&&)      = delete;
  auto operator=(const LBMBndInterface&) -> LBMBndInterface& = delete;
  auto operator=(LBMBndInterface&&) -> LBMBndInterface& = delete;

  virtual void init()                                                                                               = 0;
  virtual void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) = 0;

 private:
};

template <LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBnd_wallBB;

class LBMBnd_dummy;

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBndManager : private configuration {
 public:
  LBMBndManager()          = default;
  virtual ~LBMBndManager() = default;

  // deleted constructors not needed
  LBMBndManager(const LBMBndManager&) = delete;
  LBMBndManager(LBMBndManager&&)      = delete;
  auto operator=(const LBMBndManager&) -> LBMBndManager& = delete;
  auto operator=(LBMBndManager&&) -> LBMBndManager& = delete;

  void setupBndryCnds(const json& bndConfig, std::function<const Surface<dim(LBTYPE)>&(GString)>& bndrySurface) {
    for(const auto& [geometry, geomBndConfig] : bndConfig.items()) {
      // todo: check that geometry exists
      // todo: access only bndrySurfaces of this geometry!
      for(const auto& [surfId, surfBndConfig] : geomBndConfig.items()) {
        const auto bndType = config::required_config_value<GString>(surfBndConfig, "type");
        logger << "Adding bndCnd to surfaceId " << surfId;

        // todo: convert to switch
        if(bndType == "periodic") {
          logger << " with periodic conditions" << std::endl;
          addBndry(BndryType::Periodic, bndrySurface(surfId), surfBndConfig);
          continue;
        }
        if(bndType == "wall"){
          const GDouble tangentialV = config::opt_config_value(surfBndConfig, "tangentialVelocity", 0.0);
          if(std::abs(tangentialV) > GDoubleEps){
            logger << " wall bounce-back with tangential velocity" << std::endl;
            addBndry(BndryType::Wall_BounceBack_TangentialVelocity, bndrySurface(surfId), surfBndConfig);
          } else{
            logger << " wall bounce-back with no slip" << std::endl;
            addBndry(BndryType::Wall_BounceBack, bndrySurface(surfId), surfBndConfig);
          }
          continue;
        }
      }
    }
  }

  void addBndry(const BndryType bnd, const Surface<dim(LBTYPE)>& surf, const json& properties) {
    ASSERT(surf.size() > 0, "Invalid surface");
    switch(bnd) {
      case BndryType::Wall_BounceBack:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<LBTYPE, false>>(surf, properties));
        break;
      case BndryType::Wall_BounceBack_TangentialVelocity:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<LBTYPE, true>>(surf, properties));
        break;
      case BndryType::Periodic:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_dummy>());
        break;
      default:
        TERMM(-1, "Invalid bndry Type!");
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) {
    for(const auto& bndry : m_bndrys) {
      bndry->apply(f, fold);
    }
  }

 private:
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
  static constexpr GInt NVARS = NDIM + 1 + static_cast<GInt>(LBMethod<LBTYPE>::m_isThermal);

  std::vector<std::unique_ptr<LBMBndInterface>> m_bndrys;
};

template <LBMethodType LBTYPE>
class LBMBndCell {
 public:
  LBMBndCell() = delete;
  //  LBMBndCell() = default;
  LBMBndCell(const GInt mappedCellId, const VectorD<dim(LBTYPE)>& normal) : m_mappedCellId(mappedCellId), m_normal(normal) {}
  virtual ~LBMBndCell() = default;

  virtual void init() {
    // calculate local bndry normal
  }

  [[nodiscard]] auto mapped() const -> GInt { return m_mappedCellId; }
  [[nodiscard]] auto normal() const -> const VectorD<dim(LBTYPE)>& { return m_normal; }

  // todo: fix me
  //  deleted constructors not needed
  //  LBMBndCell(const LBMBndCell&) = delete;
  //  LBMBndCell(LBMBndCell&&)      = delete;
  //  virtual auto operator=(const LBMBndCell&) -> LBMBndCell& = delete;
  //  virtual auto operator=(LBMBndCell&&) -> LBMBndCell& = delete;

 private:
  VectorD<dim(LBTYPE)> m_normal;
  GInt                 m_mappedCellId = -1;
};

class LBMBnd_dummy : public LBMBndInterface {
 public:
  void init() override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) override {}
};

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
      if(inDirection<dim(LBTYPE)>(normal(), LBMethod<LBTYPE>::m_dirs[dist])) {
        m_bndIndex[m_noSetDists] = dist;
        ++m_noSetDists;
      }
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) {
    // iterate over the distributions that need to be set
    for(GInt id = 0; id < m_noSetDists; ++id) {
      const GInt dist         = m_bndIndex[id];
      GInt       oppositeDist = LBMethod<LBTYPE>::oppositeDist(dist);

      // standard bounceback i.e. distribution that hits the wall is reflected to the opposite distribution direction
      fold(mapped(), oppositeDist) = f(mapped(), dist);
    }

    // todo: testing
    if constexpr(TANGENTIALVELO) {
      for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
        fold(mapped(), dist) += m_tangentialVelo[dist];
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

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) override {
    for(auto& bndCell : m_bndCells) {
      bndCell.apply(f, fold);
    }
  }

 private:
  GDouble                                                m_tangentialVelo = 0.1;
  std::vector<LBMBndCell_wallBB<LBTYPE, TANGENTIALVELO>> m_bndCells;
};


#endif // LBM_LBM_BND_H
