#ifndef LBM_BND_H
#define LBM_BND_H

#include <sfcmm_common.h>
#include "common/surface.h"
#include "common/configuration.h"
#include "lbm_bnd_in_out.h"
#include "lbm_bnd_interface.h"
#include "lbm_bnd_periodic.h"
#include "lbm_bnd_wall.h"
#include "lbm_bnd_dirichlet.h"

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBnd_wallBB;

class LBMBnd_dummy;

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
class LBMBndManager : private Configuration {
 public:
  LBMBndManager()          = default;
  virtual ~LBMBndManager() = default;

  // deleted constructors not needed
  LBMBndManager(const LBMBndManager&) = delete;
  LBMBndManager(LBMBndManager&&)      = delete;
  auto operator=(const LBMBndManager&) -> LBMBndManager& = delete;
  auto operator=(LBMBndManager&&) -> LBMBndManager& = delete;

  // todo: cleanup
  void setupBndryCnds(const json& bndConfig, std::function<const Surface<DEBUG_LEVEL, dim(LBTYPE)>&(GString)>& bndrySurface) {
    for(const auto& [geometry, geomBndConfig] : bndConfig.items()) {
      // todo: check that geometry exists
      // todo: access only bndrySurfaces of this geometry!
      for(const auto& [surfId, surfBndConfig] : geomBndConfig.items()) {
        const auto bndType = config::required_config_value<GString>(surfBndConfig, "type");
        logger << "Adding bndCnd to surfaceId " << surfId;

        const Surface<DEBUG_LEVEL, dim(LBTYPE)>&              srf = bndrySurface(surfId);
        std::vector<const Surface<DEBUG_LEVEL, dim(LBTYPE)>*> bndrySrf;
        bndrySrf.emplace_back(&srf);

        // todo: convert to switch
        if(bndType == "periodic") {
          const auto surfConnected = bndrySurface(config::required_config_value<GString>(surfBndConfig, "connection"));
          bndrySrf.emplace_back(&surfConnected);
          logger << " with periodic conditions" << std::endl;
          addBndry(BndryType::Periodic, surfBndConfig, bndrySrf);
          continue;
        }
        if(bndType == "wall") {
          const GDouble tangentialV = config::opt_config_value(surfBndConfig, "tangentialVelocity", 0.0);
          if(std::abs(tangentialV) > GDoubleEps) {
            logger << " wall bounce-back with tangential velocity" << std::endl;
            addBndry(BndryType::Wall_BounceBack_TangentialVelocity, surfBndConfig, bndrySrf);
          } else {
            logger << " wall bounce-back with no slip" << std::endl;
            addBndry(BndryType::Wall_BounceBack, surfBndConfig, bndrySrf);
          }
          continue;
        }
        if(bndType == "outlet") {
          logger << " outlet bounce-back with constant pressure" << std::endl;
          addBndry(BndryType::Outlet_BounceBack_ConstPressure, surfBndConfig, bndrySrf);
          continue;
        }
        if(bndType == "inlet") {
          logger << " inlet bounce-back with constant pressure" << std::endl;
          addBndry(BndryType::Inlet_BounceBack_ConstPressure, surfBndConfig, bndrySrf);
          continue;
        }
        if(bndType == "dirichlet") {
          logger << " dirichlet boundary condition using neem" << std::endl;
          addBndry(BndryType::Dirichlet_NEEM, surfBndConfig, bndrySrf);
          continue;
        }
        TERMM(-1, "Invalid bndCndType: " + bndType);
      }
    }
  }

  // todo: merge with the function above
  void addBndry(const BndryType bnd, const json& properties, const std::vector<const Surface<DEBUG_LEVEL, dim(LBTYPE)>*>& surf) {
    ASSERT(surf[0]->size() > 0, "Invalid surface");

    // actually generate a boundary or just create a dummy e.g. the boundary is handled in a non standard way!
    const auto generateBndry = config::opt_config_value<GBool>(properties, "generateBndry", true);

    if(generateBndry) {
      switch(bnd) {
        case BndryType::Wall_BounceBack:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<DEBUG_LEVEL, LBTYPE, false>>(surf[0], properties));
          break;
        case BndryType::Wall_BounceBack_TangentialVelocity:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<DEBUG_LEVEL, LBTYPE, true>>(surf[0], properties));
          break;
        case BndryType::Inlet_BounceBack_ConstPressure:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_InOutBB<DEBUG_LEVEL, LBTYPE>>(surf[0], properties));
          TERMM(-1, "Broken");
          break;
        case BndryType::Outlet_BounceBack_ConstPressure:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_InOutBB<DEBUG_LEVEL, LBTYPE>>(surf[0], properties));
          TERMM(-1, "Broken");
          break;
        case BndryType::Periodic:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_Periodic<DEBUG_LEVEL, LBTYPE>>(surf[0], surf[1], properties));
          break;
        case BndryType::Dirichlet_NEEM:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_DirichletNEEM<DEBUG_LEVEL, LBTYPE, EQ>>(surf[0], properties));
          break;
        default:
          TERMM(-1, "Invalid bndry Type!");
      }
    } else {
      logger << "Generated dummy boundary!" << std::endl;
      m_bndrys.emplace_back(std::make_unique<LBMBnd_dummy>());
    }
  }

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                const std::function<GDouble&(GInt, GInt)>& vars) {
    for(const auto& bndry : m_bndrys) {
      bndry->preApply(f, fold, vars);
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& vars) {
    for(const auto& bndry : m_bndrys) {
      bndry->apply(f, fold, vars);
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
  LBMBnd_dummy()           = default;
  ~LBMBnd_dummy() override = default;

  //  void init() override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}
  void apply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
             const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}
};
#endif // LBM_BND_H
