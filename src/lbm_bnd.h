#ifndef LBM_BND_H
#define LBM_BND_H

#include <sfcmm_common.h>
#include "common/surface.h"
#include "configuration.h"
#include "lbm_bnd_in_out.h"
#include "lbm_bnd_interface.h"
#include "lbm_bnd_wall.h"

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
        if(bndType == "wall") {
          const GDouble tangentialV = config::opt_config_value(surfBndConfig, "tangentialVelocity", 0.0);
          if(std::abs(tangentialV) > GDoubleEps) {
            logger << " wall bounce-back with tangential velocity" << std::endl;
            addBndry(BndryType::Wall_BounceBack_TangentialVelocity, bndrySurface(surfId), surfBndConfig);
          } else {
            logger << " wall bounce-back with no slip" << std::endl;
            addBndry(BndryType::Wall_BounceBack, bndrySurface(surfId), surfBndConfig);
          }
          continue;
        }
        if(bndType == "outlet"){
          logger << " outlet bounce-back with constant pressure" << std::endl;
          addBndry(BndryType::Outlet_BounceBack_ConstPressure, bndrySurface(surfId), surfBndConfig);
          continue;
        }
        if(bndType == "inlet"){
          logger << " inlet bounce-back with constant pressure" << std::endl;
          addBndry(BndryType::Inlet_BounceBack_ConstPressure, bndrySurface(surfId), surfBndConfig);
          continue;
        }
        TERMM(-1, "Invalid bndCndType: " + bndType);
      }
    }
  }

  //todo: merge with the function above
  void addBndry(const BndryType bnd, const Surface<dim(LBTYPE)>& surf, const json& properties) {
    ASSERT(surf.size() > 0, "Invalid surface");
    switch(bnd) {
      case BndryType::Wall_BounceBack:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<LBTYPE, false>>(surf, properties));
        break;
      case BndryType::Wall_BounceBack_TangentialVelocity:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<LBTYPE, true>>(surf, properties));
        break;
      case BndryType::Inlet_BounceBack_ConstPressure:
      case BndryType::Outlet_BounceBack_ConstPressure:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_InOutBB<LBTYPE>>(surf, properties));
        break;
      case BndryType::Periodic:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_dummy>());
        break;
      default:
        TERMM(-1, "Invalid bndry Type!");
    }
  }

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) {
    for(const auto& bndry : m_bndrys) {
      bndry->preApply(f, fold);
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

class LBMBnd_dummy : public LBMBndInterface {
 public:
  void init() override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) override {}
  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) override {}
};
#endif // LBM_BND_H
