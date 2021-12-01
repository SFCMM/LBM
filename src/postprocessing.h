#ifndef LBM_POSTPROCESSING_H
#define LBM_POSTPROCESSING_H

#include <lbm/lbm_postprocessing.h>
#include <utility>
#include "configuration.h"
#include "postprocessing_func.h"

namespace pp {
enum class HOOK { ATSTART, BEFORETIMESTEP, AFTERTIMESTEP, ATEND, NUM };

inline auto getHook(const GString& hookStr) -> HOOK {
  if(hookStr == "atStart") {
    return HOOK::ATSTART;
  }
  if(hookStr == "beforeTimestep") {
    return HOOK::BEFORETIMESTEP;
  }
  if(hookStr == "afterTimestep") {
    return HOOK::AFTERTIMESTEP;
  }
  if(hookStr == "atEnd") {
    return HOOK::ATEND;
  }
  TERMM(-1, "Invalid hook");
}

inline auto getHook(const HOOK hook) -> GString {
  if(hook == HOOK::ATSTART) {
    return "atStart";
  }
  if(hook == HOOK::BEFORETIMESTEP) {
    return "beforeTimestep";
  }
  if(hook == HOOK::AFTERTIMESTEP) {
    return "afterTimestep";
  }
  if(hook == HOOK::ATEND) {
    return "atEnd";
  }
  TERMM(-1, "Invalid hook");
}
} // namespace pp

template <Debug_Level DEBUG_LEVEL, GInt NDIM, SolverType SOLVT, class SOLV>
class Postprocess {
 public:
  Postprocess() = default;

  void init() {
    if(!m_postprocessing){
      //nothing to do not active
      return;
    }

    if(SOLVT == SolverType::LBM) {
      for(const auto& obj : m_conf->getAllObjects()) {
        const auto         typeId = m_conf->required_config_value<GString>(obj, "type");
        const auto         exeAt  = m_conf->required_config_value<GString>(obj, "execute");
        const pp::FuncType ppType = pp::getFuncType(typeId);
        const GInt         atHook = static_cast<GInt>(pp::getHook(exeAt));

        logger << "PP: Adding postprocessing function " << typeId << " for LBM solver at hook" << exeAt << std::endl;

        switch(ppType) {
          case pp::FuncType::LINE:
            m_postProcessFunc[atHook].emplace_back(
                std::make_unique<LBMPostprocessFunctionLine<NDIM>>(m_conf->getObject(obj), solver()->getCartesianGridData()));
            break;
          default:
            TERMM(-1, "Invalid postprocessing function type");
        }
      }
    } else {
      TERMM(-1, "Invalid solver type");
    }
  }

  void setConfAccessor(std::shared_ptr<ConfigurationAccess> ppConf) {
    m_conf = std::move(ppConf);
    if(m_conf != nullptr) {
      logger << "Postprocessing is ACTIVE" << std::endl;
      m_postprocessing = true;
    }
  }

  void executePostprocess(const pp::HOOK hook) {
    if(m_postprocessing) {
      const GInt ppHookId = static_cast<GInt>(hook);
      if(!m_postProcessFunc.at(ppHookId).empty()) {
        logger << "Executing postprocessing at hook:" << pp::getHook(hook) << std::endl;
        cerr0 << "Executing postprocessing at hook:" << pp::getHook(hook) << std::endl;
      }
    }
  }

 private:


  inline auto solver() -> SOLV*{
    return static_cast<SOLV*>(this);
  }

  std::array<std::vector<std::unique_ptr<PostprocessFunctionInterface<NDIM>>>, static_cast<GInt>(pp::HOOK::NUM)> m_postProcessFunc;
  GBool                                                                                                          m_postprocessing = false;
  std::shared_ptr<ConfigurationAccess>                                                                           m_conf           = nullptr;
};
#endif // LBM_POSTPROCESSING_H
