#ifndef LBM_BND_INTERFACE_H
#define LBM_BND_INTERFACE_H
#include "lbm/constants.h"
#include "sfcmm_common.h"

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

class LBMBnd_dummy : public LBMBndInterface {
 public:
  LBMBnd_dummy()           = default;
  ~LBMBnd_dummy() override = default;

  // deleted constructors not needed
  LBMBnd_dummy(const LBMBnd_dummy&)                    = delete;
  LBMBnd_dummy(LBMBnd_dummy&&)                         = delete;
  auto operator=(const LBMBnd_dummy&) -> LBMBnd_dummy& = delete;
  auto operator=(LBMBnd_dummy&&) -> LBMBnd_dummy&      = delete;

  void initCnd(const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}
  void apply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
             const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}
};
#endif // LBM_BND_INTERFACE_H
