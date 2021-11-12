#ifndef LBM_LBM_BND_IN_OUT_H
#define LBM_LBM_BND_IN_OUT_H
#include "lbm_bnd_interface.h"
#include "lbm_constants.h"

template <LBMethodType LBTYPE>
class LBMBnd_InOutBB : public LBMBndInterface {
 public:
  LBMBnd_InOutBB(const Surface<dim(LBTYPE)>& surf, const json& properties) {
    for(const GInt cellId : surf.getCellList()) {
      m_bndCells.emplace_back(cellId);
    }
    LBMBnd_InOutBB<LBTYPE>::init();
  }
  ~LBMBnd_InOutBB() override = default;

  // deleted constructors not needed
  LBMBnd_InOutBB(const LBMBnd_InOutBB&) = delete;
  LBMBnd_InOutBB(LBMBnd_InOutBB&&)      = delete;
  auto operator=(const LBMBnd_InOutBB&) -> LBMBnd_InOutBB& = delete;
  auto operator=(LBMBnd_InOutBB&&) -> LBMBnd_InOutBB& = delete;

  void init() override {

  }

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) override {

  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) override {

  }



 private:
  GDouble           m_pressure = 1.0;
  std::vector<GInt> m_bndCells;
};

#endif // LBM_LBM_BND_IN_OUT_H
