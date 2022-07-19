#ifndef LBM_BND_INTERFACE_H
#define LBM_BND_INTERFACE_H
#include <sfcmm_common.h>
#include "constants.h"

class LBMBndInterface {
 public:
  LBMBndInterface()          = default;
  virtual ~LBMBndInterface() = default;

  // deleted constructors not needed
  LBMBndInterface(const LBMBndInterface&)                    = delete;
  LBMBndInterface(LBMBndInterface&&)                         = delete;
  auto operator=(const LBMBndInterface&) -> LBMBndInterface& = delete;
  auto operator=(LBMBndInterface&&) -> LBMBndInterface&      = delete;

  virtual void initCnd(const std::function<GDouble&(GInt, GInt)>& vars)                                                  = 0;
  virtual void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                        const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) = 0;
  virtual void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                     const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars)    = 0;

 private:
};
#endif // LBM_BND_INTERFACE_H
