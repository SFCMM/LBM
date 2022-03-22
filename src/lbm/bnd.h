#ifndef LBM_BND_H
#define LBM_BND_H

#include <sfcmm_common.h>
#include "bnd_dirichlet.h"
#include "bnd_neumann.h"
#include "bnd_in_out.h"
#include "bnd_interface.h"
#include "bnd_periodic.h"
#include "bnd_wall.h"
#include "common/configuration.h"
#include "common/surface.h"

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, GBool TANGENTIALVELO>
class LBMBnd_wallBB;

class LBMBnd_dummy;

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ>
class LBMBndManager : private Configuration {
 private:
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
  static constexpr GInt NVARS = noVars<LBTYPE>(EQ);

 public:
  LBMBndManager()          = default;
  virtual ~LBMBndManager() = default;

  // deleted constructors not needed
  LBMBndManager(const LBMBndManager&) = delete;
  LBMBndManager(LBMBndManager&&)      = delete;
  auto operator=(const LBMBndManager&) -> LBMBndManager& = delete;
  auto operator=(LBMBndManager&&) -> LBMBndManager& = delete;


  /// Setup the boundary conditions.
  /// \param bndConfig Configuration of the boundaries
  /// \param bndrySurface List of the boundaries
  void setupBndryCnds(const json& bndConfig, std::function<Surface<DEBUG_LEVEL, dim(LBTYPE)>&(GString)>& bndrySurface) {
    for(const auto& [geometry, geomBndConfig] : bndConfig.items()) {
      // todo: cleanup
      // todo: check that geometry exists
      // todo: access only bndrySurfaces of this geometry!
      const GInt noBnds = geomBndConfig.size();
      for(const auto& [surfId, surfBndConfig] : geomBndConfig.items()) {
        GString surfIdName = surfId;
        if(noBnds > 1) {
          surfIdName = geometry + "_" + surfIdName;
        }
        const auto bndType = config::required_config_value<GString>(surfBndConfig, "type");
        logger << "Adding bndCnd to surfaceId " << surfIdName;

        Surface<DEBUG_LEVEL, dim(LBTYPE)>& srf = bndrySurface(surfIdName);
        if(srf.size() == 0) {
          continue;
        }
        std::vector<Surface<DEBUG_LEVEL, dim(LBTYPE)>*> bndrySrf;
        bndrySrf.emplace_back(&srf);

        // todo: convert to switch
        if(bndType == "periodic") {
          auto surfConnected = bndrySurface(config::required_config_value<GString>(surfBndConfig, "connection"));
          ASSERT(surfConnected.size() > 0, "Invalid surface");
          bndrySrf.emplace_back(&surfConnected);
          logger << " with periodic conditions" << std::endl;
          addBndry(BndryType::Periodic, surfBndConfig, bndrySrf);
          continue;
        }
        if(bndType == "wall") {
          const auto model = config::required_config_value<GString>(surfBndConfig, "model");
          if(model == "bounceback") {
            const GDouble tangentialV = config::opt_config_value(surfBndConfig, "tangentialVelocity", 0.0);
            if(std::abs(tangentialV) > GDoubleEps) {
              logger << " wall bounce-back with tangential velocity" << std::endl;
              addBndry(BndryType::Wall_BounceBack_TangentialVelocity, surfBndConfig, bndrySrf);
            } else {
              logger << " wall bounce-back with no slip" << std::endl;
              addBndry(BndryType::Wall_BounceBack, surfBndConfig, bndrySrf);
            }
          } else if(model == "equilibrium") {
            logger << " wall equilibrium wet node model" << std::endl;
            addBndry(BndryType::Wall_Equilibrium, surfBndConfig, bndrySrf);
          } else if(model == "neem") {
            logger << " wall nonequilibrium extrapolation wet node model" << std::endl;
            addBndry(BndryType::Wall_NEEM, surfBndConfig, bndrySrf);
          } else if(model == "nebb") {
            logger << " wall nonequilibrium bounce back wet node model" << std::endl;
            addBndry(BndryType::Wall_NEBB, surfBndConfig, bndrySrf);
          } else {
            TERMM(-1, "Invalid wall boundary model: " + model);
          }
          continue;
        }
        if(bndType == "pressure") {
          logger << " open constant pressure anti bounce-back condition" << std::endl;
          addBndry(BndryType::Open_AntiBounceBack_ConstPressure, surfBndConfig, bndrySrf);
          continue;
        }
        if(bndType == "outlet") {
          logger << " outlet anti bounce-back with constant pressure" << std::endl;
          addBndry(BndryType::Outlet_AntiBounceBack_ConstPressure, surfBndConfig, bndrySrf);
          continue;
        }
        if(bndType == "inlet") {
          logger << " inlet anti bounce-back with constant pressure" << std::endl;
          addBndry(BndryType::Inlet_AntiBounceBack_ConstPressure, surfBndConfig, bndrySrf);
          continue;
        }
        if(bndType == "dirichlet") {
          const auto model = config::required_config_value<GString>(surfBndConfig, "model");
          if(model == "bounceback") {
            logger << " dirichlet boundary condition using bb" << std::endl;
            addBndry(BndryType::Dirichlet_BounceBack, surfBndConfig, bndrySrf);
          } else if(model == "neem") {
            logger << " dirichlet boundary condition using neem" << std::endl;
            addBndry(BndryType::Dirichlet_NEEM, surfBndConfig, bndrySrf);
          }
          continue;
        }
        if(bndType == "neumann") {
          const auto model = config::required_config_value<GString>(surfBndConfig, "model");
          //          if(model == "bounceback") {
          //            logger << " neumann boundary condition using bb" << std::endl;
          //            addBndry(BndryType::Neumann_BounceBack, surfBndConfig, bndrySrf);
          //          } else
          if(model == "neem") {
            logger << " neumann boundary condition using neem" << std::endl;
            addBndry(BndryType::Neumann_NEEM, surfBndConfig, bndrySrf);
          }
          continue;
        }
        TERMM(-1, "Invalid bndCndType: " + bndType);
      }
    }
  }

  void initCndBnd(const std::function<GDouble&(GInt, GInt)>& vars) {
    for(const auto& bndry : m_bndrys) {
      bndry->initCnd(vars);
    }
  }

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) {
    for(const auto& bndry : m_bndrys) {
      bndry->preApply(f, fold, feq, vars);
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) {
    for(const auto& bndry : m_bndrys) {
      bndry->apply(f, fold, feq, vars);
    }
  }

 private:
  // todo: merge with the function above
  void addBndry(const BndryType bnd, const json& properties, const std::vector<Surface<DEBUG_LEVEL, dim(LBTYPE)>*>& surf) {
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
        case BndryType::Wall_Equilibrium:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_wallEq<DEBUG_LEVEL, LBTYPE>>(surf[0], properties));
          break;
        case BndryType::Wall_NEEM:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_wallNEEM<DEBUG_LEVEL, LBTYPE, EQ>>(surf[0], properties));
          break;
        case BndryType::Wall_NEBB:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_wallNEBB<DEBUG_LEVEL, LBTYPE>>(surf[0], properties));
          break;
        case BndryType::Open_AntiBounceBack_ConstPressure:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_Pressure<DEBUG_LEVEL, LBTYPE>>(surf[0], properties));
          break;
        case BndryType::Inlet_AntiBounceBack_ConstPressure:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_Pressure<DEBUG_LEVEL, LBTYPE>>(surf[0], properties));
          TERMM(-1, "Broken");
          //          break;
        case BndryType::Outlet_AntiBounceBack_ConstPressure:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_Pressure<DEBUG_LEVEL, LBTYPE>>(surf[0], properties));
          TERMM(-1, "Broken");
          //          break;
        case BndryType::Periodic:
          ASSERT(surf[1]->size() > 0, "Invalid connected surface");
          surf[0]->setProperty(CellProperties::periodic, true);
          surf[1]->setProperty(CellProperties::periodic, true);
          m_bndrys.emplace_back(std::make_unique<LBMBnd_Periodic<DEBUG_LEVEL, LBTYPE>>(surf[0], surf[1], properties));
          break;
        case BndryType::Dirichlet_NEEM:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_DirichletNEEM<DEBUG_LEVEL, LBTYPE, EQ>>(surf[0], properties));
          break;
        case BndryType::Dirichlet_BounceBack:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_DirichletBB<DEBUG_LEVEL, LBTYPE, EQ>>(surf[0], properties));
          break;
        case BndryType::Neumann_NEEM:
          m_bndrys.emplace_back(std::make_unique<LBMBnd_NeumannNEEM<DEBUG_LEVEL, LBTYPE, EQ>>(surf[0], properties));
          break;
          // todo:implement
          //        case BndryType::Dirichlet_BounceBack:
          //          m_bndrys.emplace_back(std::make_unique<LBMBnd_DirichletBB<DEBUG_LEVEL, LBTYPE, EQ>>(surf[0], properties));
          //          break;
        default:
          TERMM(-1, "Invalid bndry Type!");
      }
    } else {
      logger << "Generated dummy boundary!" << std::endl;
      m_bndrys.emplace_back(std::make_unique<LBMBnd_dummy>());
    }
  }

  std::vector<std::unique_ptr<LBMBndInterface>> m_bndrys;
};

class LBMBnd_dummy : public LBMBndInterface {
 public:
  LBMBnd_dummy()           = default;
  ~LBMBnd_dummy() override = default;

  // deleted constructors not needed
  LBMBnd_dummy(const LBMBnd_dummy&) = delete;
  LBMBnd_dummy(LBMBnd_dummy&&)      = delete;
  auto operator=(const LBMBnd_dummy&) -> LBMBnd_dummy& = delete;
  auto operator=(LBMBnd_dummy&&) -> LBMBnd_dummy& = delete;

  void initCnd(const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}
  void apply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
             const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}
};
#endif // LBM_BND_H
