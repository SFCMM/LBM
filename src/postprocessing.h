#ifndef LBM_POSTPROCESSING_H
#define LBM_POSTPROCESSING_H

#include <utility>
#include "configuration.h"
#include "postprocessing_func.h"

namespace pp {
enum class HOOK { ATSTART, BEFORETIMESTEP, AFTERTIMESTEP, ATEND, NUM };

inline auto post(const GString& hookStr) -> HOOK {
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
}

inline auto post(const HOOK hook) -> GString {
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

template <Debug_Level DEBUG_LEVEL, GInt NDIM, SolverType SOLVT>
class Postprocess {
 public:
  Postprocess() = default;

  void setConfAccessor(ConfigurationAccess* ppConf) {
    m_conf = ppConf;
    if(m_conf != nullptr) {
      m_postprocessing = true;
      init();
    }
  }

  void executePostprocess(const pp::HOOK hook) {
    if(m_postprocessing) {
      const GInt ppHookId = static_cast<GInt>(hook);
      if(!m_postProcessFunc.at(ppHookId).empty()) {
        logger << "Executing postprocessing at hook:" << pp::post(hook) << std::endl;
        cerr0 << "Executing postprocessing at hook:" << pp::post(hook) << std::endl;
      }
    }
  }

 private:
  void init() {
    if(SOLVT == SolverType::LBM){

    } else{
      TERMM(-1, "Invalid solver type");
    }
  }

  std::array<std::vector<std::unique_ptr<PostprocessFunctionExecutorInterface<NDIM>>>, static_cast<GInt>(pp::HOOK::NUM)> m_postProcessFunc;
  GBool                m_postprocessing = false;
  ConfigurationAccess* m_conf           = nullptr;
};
#endif // LBM_POSTPROCESSING_H
