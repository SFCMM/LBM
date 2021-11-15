#ifndef LBM_LBM_BND_IN_OUT_H
#define LBM_LBM_BND_IN_OUT_H
#include "lbm_bnd_interface.h"
#include "lbm_constants.h"
#include "pv.h"

template <LBMethodType LBTYPE>
class LBMBnd_InOutBB : public LBMBndInterface {
 public:
  LBMBnd_InOutBB(const Surface<dim(LBTYPE)>& surf, const json& properties) {
    for(const GInt cellId : surf.getCellList()) {
      m_bndCells.emplace_back(cellId);
    }
    m_normal = surf.normal(0);
    LBMBnd_InOutBB<LBTYPE>::init();
  }
  ~LBMBnd_InOutBB() override = default;

  // deleted constructors not needed
  LBMBnd_InOutBB(const LBMBnd_InOutBB&) = delete;
  LBMBnd_InOutBB(LBMBnd_InOutBB&&)      = delete;
  auto operator=(const LBMBnd_InOutBB&) -> LBMBnd_InOutBB& = delete;
  auto operator=(LBMBnd_InOutBB&&) -> LBMBnd_InOutBB& = delete;

  void init() override {}

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& vars) override {
    const GDouble targetRho = 1.075;
    for(const auto cellId : m_bndCells) {
      const GDouble rho = vars(cellId, PV::rho<dim(LBTYPE)>());
      for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
        //        const GDouble oldV = fold(cellId, dist);
        //        const GDouble rho  = 0.001;
        //        fold(cellId, dist) = f(cellId, dist) + LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho); //slightly low
        if(m_normal[0] > 0) {
          // outlet
          //          fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (1.0 - rho);
          //          if(dist == 4 || dist == 1 || dist == 5 || dist == 8) {
          //            fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
          //          } else
          //              if(dist == 7 || dist == 0 || dist == 6) {
          //            fold(cellId, dist) -= LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
          //          }
          //          if(dist == 8){
          fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (1.0 - rho);
          //          }
        } else {
          if(cellId == 341 || cellId == 0){
            fold(cellId, 8) += LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
            continue;
          }
          // inlet
          if(dist == 4 || dist == 1 || dist == 5 || dist == 8) {
            //          if(dist == 8) {
//            if((dist == 4) && cellId == 341){
//              continue;
//            }
            fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
            //            fold(cellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
            //          }
          } else if(dist == 0) {
            // outside direction
            //            fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
            fold(cellId, dist) =
                -fold(cellId, LBMethod<LBTYPE>::oppositeDist(dist)) + LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
          } else if(dist == 7 || dist == 6 || dist == 3 || dist == 2) {
//            if((dist == 3 || dist==7) && cellId == 341){
//              continue;
//            }
            fold(cellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
          }
        }
      }
    }
    //    TERMM(-1, "TEsts");
  }


 private:
  GDouble              m_pressure = 1.0;
  VectorD<dim(LBTYPE)> m_normal;
  std::vector<GInt>    m_bndCells;
};

#endif // LBM_LBM_BND_IN_OUT_H
