#ifndef LBM_BND_H
#define LBM_BND_H

#include "bnd_dirichlet.h"
#include "bnd_in_out.h"
#include "bnd_interface.h"
#include "bnd_neumann.h"
#include "bnd_periodic.h"
#include "bnd_wall.h"
#include "common/configuration.h"
#include "common/surface.h"
#include "sfcmm_common.h"

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
  LBMBndManager(const LBMBndManager&)                    = delete;
  LBMBndManager(LBMBndManager&&)                         = delete;
  auto operator=(const LBMBndManager&) -> LBMBndManager& = delete;
  auto operator=(LBMBndManager&&) -> LBMBndManager&      = delete;

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


  /// Setup the boundary conditions.
  /// \param bndConfig Configuration of the boundaries
  /// \param bndrySurface List of the boundaries
  void setupBndryCnds(const json& bndConfig, std::function<Surface<DEBUG_LEVEL, dim(LBTYPE)>&(GString)>& bndrySurface) {
    // iterate over all geometries and the associated boundary setup for this geometry
    for(const auto& [geometry, geomBndConfig] : bndConfig.items()) {
      // todo: check that geometry exists
      const GInt noBnds = static_cast<GInt>(geomBndConfig.size());

      // iterate over all boundary setup configuration items which define the boundary condition to be used for each of the surfaces
      for(const auto& [surfId, surfBndConfig] : geomBndConfig.items()) {
        // give each surface a *unique* name, e.g. plane_x+ the plane's surfaces in x+ direction
        const GString surfIdName = noBnds > 1 ? geometry + "_" + surfId : geometry;

        Surface<DEBUG_LEVEL, dim(LBTYPE)>& srf = bndrySurface(surfIdName);
        if(srf.size() == 0) {
          logger << "WARNING: Skipping " << surfIdName << " no valid cells!" << std::endl;
          continue;
        }


        logger << "Adding bndCnd to surfaceId " << surfIdName;
        const auto bndType = config::required_config_value<GString>(surfBndConfig, "type");
        // todo: convert to switch
        if(bndType == "periodic") {
          auto surfConnected = bndrySurface(config::required_config_value<GString>(surfBndConfig, "connection"));
          addPeriodicBndry(BndryType::Periodic, surfBndConfig, srf, surfConnected);
          continue;
        }
        if(bndType == "wall") {
          addWallBndry(srf, surfBndConfig);
          continue;
        }
        if(bndType == "pressure") {
          logger << " open constant pressure anti bounce-back condition" << std::endl;
          addBndry(BndryType::Open_AntiBounceBack_ConstPressure, surfBndConfig, srf);
          continue;
        }
        if(bndType == "outlet") {
          logger << " outlet anti bounce-back with constant pressure" << std::endl;
          addBndry(BndryType::Outlet_AntiBounceBack_ConstPressure, surfBndConfig, srf);
          continue;
        }
        if(bndType == "inlet") {
          logger << " inlet anti bounce-back with constant pressure" << std::endl;
          addBndry(BndryType::Inlet_AntiBounceBack_ConstPressure, surfBndConfig, srf);
          continue;
        }
        if(bndType == "dirichlet") {
          const auto model = config::required_config_value<GString>(surfBndConfig, "model");
          if(model == "bounceback") {
            logger << " dirichlet boundary condition using bb" << std::endl;
            addBndry(BndryType::Dirichlet_BounceBack, surfBndConfig, srf);
          } else if(model == "neem") {
            logger << " dirichlet boundary condition using neem" << std::endl;
            addBndry(BndryType::Dirichlet_NEEM, surfBndConfig, srf);
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
            addBndry(BndryType::Neumann_NEEM, surfBndConfig, srf);
          }
          continue;
        }
        TERMM(-1, "Invalid bndCndType: " + bndType);
      }
    }
  }


 private:
  auto handleDummyBnd(Surface<DEBUG_LEVEL, dim(LBTYPE)>& /*surf*/, const json& properties) -> GBool {
    // actually generate a boundary or just create a dummy e.g. the boundary is handled in a non standard way!
    const auto generateBndry = config::opt_config_value<GBool>(properties, "generateBndry", true);

    if(!generateBndry) {
      logger << "Generated dummy boundary! (generateBndry = false)" << std::endl;
      m_bndrys.emplace_back(std::make_unique<LBMBnd_dummy>());
      return true;
    }
    return false;
  }

  void addBndry(const BndryType bnd, const json& properties, Surface<DEBUG_LEVEL, dim(LBTYPE)>& surf) {
    if(handleDummyBnd(surf, properties)) {
      return;
    }

    switch(bnd) {
      case BndryType::Open_AntiBounceBack_ConstPressure:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_Pressure<DEBUG_LEVEL, LBTYPE>>(&surf, properties));
        break;
      case BndryType::Inlet_AntiBounceBack_ConstPressure:
        // todo: fix
        m_bndrys.emplace_back(std::make_unique<LBMBnd_Pressure<DEBUG_LEVEL, LBTYPE>>(&surf, properties));
        TERMM(-1, "Broken");
        //          break;
      case BndryType::Outlet_AntiBounceBack_ConstPressure:
        // todo: fix
        m_bndrys.emplace_back(std::make_unique<LBMBnd_Pressure<DEBUG_LEVEL, LBTYPE>>(&surf, properties));
        TERMM(-1, "Broken");
        //          break;
      case BndryType::Dirichlet_NEEM:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_DirichletNEEM<DEBUG_LEVEL, LBTYPE, EQ>>(&surf, properties));
        break;
      case BndryType::Dirichlet_BounceBack:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_DirichletBB<DEBUG_LEVEL, LBTYPE, EQ>>(&surf, properties));
        break;
      case BndryType::Neumann_NEEM:
        m_bndrys.emplace_back(std::make_unique<LBMBnd_NeumannNEEM<DEBUG_LEVEL, LBTYPE, EQ>>(&surf, properties));
        break;
        // todo:implement
        //        case BndryType::Dirichlet_BounceBack:
        //          m_bndrys.emplace_back(std::make_unique<LBMBnd_DirichletBB<DEBUG_LEVEL, LBTYPE, EQ>>(surf[0], properties));
        //          break;
      default:
        TERMM(-1, "Invalid bndry Type!");
    }
  }

  void addPeriodicBndry(const BndryType bnd, const json& properties, Surface<DEBUG_LEVEL, dim(LBTYPE)>& surfA,
                        Surface<DEBUG_LEVEL, dim(LBTYPE)>& surfB) {
    if(handleDummyBnd(surfA, properties)) {
      return;
    }
    logger << " with periodic conditions" << std::endl;

    ASSERT(surfA.size() > 0, "Invalid surface");
    ASSERT(surfB.size() > 0, "Invalid surface");


    switch(bnd) {
      case BndryType::Periodic:
        surfA.setProperty(CellProperties::periodic, true);
        surfB.setProperty(CellProperties::periodic, true);
        m_bndrys.emplace_back(std::make_unique<LBMBnd_Periodic<DEBUG_LEVEL, LBTYPE>>(&surfA, &surfB, properties));
        break;
      default:
        TERMM(-1, "Invalid bndry Type!");
    }
  }


  void addWallBndry(Surface<DEBUG_LEVEL, dim(LBTYPE)>& surf, const json& properties) {
    if(handleDummyBnd(surf, properties)) {
      return;
    }
    const auto model = config::required_config_value<GString>(properties, "model");

    // todo: convert to switch
    if(model == "bounceback") {
      const GDouble tangentialV = config::opt_config_value(properties, "tangentialVelocity", 0.0);
      if(std::abs(tangentialV) > GDoubleEps) {
        logger << " wall bounce-back with tangential velocity" << std::endl;
        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<DEBUG_LEVEL, LBTYPE, true>>(&surf, properties));
      } else {
        logger << " wall bounce-back with no slip" << std::endl;
        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<DEBUG_LEVEL, LBTYPE, false>>(&surf, properties));
      }
    } else if(model == "equilibrium") {
      logger << " wall equilibrium wet node model" << std::endl;
      m_bndrys.emplace_back(std::make_unique<LBMBnd_wallEq<DEBUG_LEVEL, LBTYPE>>(&surf, properties));
    } else if(model == "neem") {
      logger << " wall nonequilibrium extrapolation wet node model" << std::endl;
      m_bndrys.emplace_back(std::make_unique<LBMBnd_wallNEEM<DEBUG_LEVEL, LBTYPE, EQ>>(&surf, properties));
    } else if(model == "nebb") {
      logger << " wall nonequilibrium bounce back wet node model" << std::endl;
      m_bndrys.emplace_back(std::make_unique<LBMBnd_wallNEBB<DEBUG_LEVEL, LBTYPE>>(&surf, properties));
    } else {
      TERMM(-1, "Invalid wall boundary model: " + model);
    }

    //    switch(bnd) {
    //      case BndryType::Wall_BounceBack:
    //        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<DEBUG_LEVEL, LBTYPE, false>>(&surf, properties));
    //        break;
    //      case BndryType::Wall_BounceBack_TangentialVelocity:
    //        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallBB<DEBUG_LEVEL, LBTYPE, true>>(&surf, properties));
    //        break;
    //      case BndryType::Wall_Equilibrium:
    //        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallEq<DEBUG_LEVEL, LBTYPE>>(&surf, properties));
    //        break;
    //      case BndryType::Wall_NEEM:
    //        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallNEEM<DEBUG_LEVEL, LBTYPE, EQ>>(&surf, properties));
    //        break;
    //      case BndryType::Wall_NEBB:
    //        m_bndrys.emplace_back(std::make_unique<LBMBnd_wallNEBB<DEBUG_LEVEL, LBTYPE>>(&surf, properties));
    //        break;
    //      default:
    //        TERMM(-1, "Invalid bndry Type!");
    //    }
  }

  std::vector<std::unique_ptr<LBMBndInterface>> m_bndrys;
};


#endif // LBM_BND_H
