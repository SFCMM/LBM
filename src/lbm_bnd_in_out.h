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

  void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                const std::function<GDouble&(GInt, GInt)>& vars) override {
    static constexpr std::array<GDouble, 9> cx = {-1, 1, 0, 0, 1, 1, -1, -1, 0};
    static constexpr std::array<GDouble, 9> cy = {0, 0, -1, 1, 1, -1, -1, 1, 0};

    //    const GDouble diffRho   = 0.15;
    const GDouble targetRho = 1.008;
    for(const auto cellId : m_bndCells) {
      const GDouble rho = vars(cellId, PV::rho<dim(LBTYPE)>());
      for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
        if(m_normal[0] < 0) {
          if(cellId == 341 || cellId == 0) {
            f(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
            fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
            continue;
          }

          // inlet
          f(cellId, dist) =
              LBMethod<LBTYPE>::m_weights[dist] * targetRho * (1 + 3 * (vars(cellId, PV::velocitiy(0)) * cx[dist] + 0 * cy[dist]));
          fold(cellId, dist) =
              LBMethod<LBTYPE>::m_weights[dist] * targetRho * (1 + 3 * (vars(cellId, PV::velocitiy(0)) * cx[dist] + 0 * cy[dist]));


          //          f(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * diffRho;
        } else {
//          if(cellId == 682 || cellId == 1023) {
//            f(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (0.9951 - rho);
//            fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (0.9951 - rho);
//            continue;
//          }
//          f(cellId, dist)    = LBMethod<LBTYPE>::m_weights[dist] * (0.9951 + 3 * (vars(cellId, PV::velocitiy(0)) * cx[dist] + 0 *
//                                                                                                                               cy[dist]));
//          fold(cellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * (0.9951 + 3 * (vars(cellId, PV::velocitiy(0)) * cx[dist] + 0 *
//                                                                                                                                  cy[dist]));
        }
      }
    }
  }

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& vars) override {
    static constexpr std::array<GDouble, 9> cx = {-1, 1, 0, 0, 1, 1, -1, -1, 0};
    static constexpr std::array<GDouble, 9> cy = {0, 0, -1, 1, 1, -1, -1, 1, 0};
    //    const GDouble targetRho = 1.075;
    //todo: both rho and velocity need to be recalculated move to a separate function
    //todo: function to recalculate feq
    for(const auto cellId : m_bndCells) {
      GDouble u = 0;
      GDouble rho = 0;
      for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
        u += fold(cellId, dist) * cx[dist];
        rho += fold(cellId, dist);
      }
      if(m_normal[0] > 0) {
        for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
          // set zero y-velocity
          fold(cellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * (1.0 + 3 * (u * cx[dist] + 0 * cy[dist]));
        }
      } else {
        for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
          // set zero y-velocity
          fold(cellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * rho * (1.0 + 3 * (u * cx[dist] + 0 * cy[dist]));
        }
      }
    }
    //      const GDouble rho = vars(cellId, PV::rho<dim(LBTYPE)>());
    //      for(GInt dist = 0; dist < noDists(LBTYPE); ++dist) {
    //        const GDouble oldV = fold(cellId, dist);
    //        const GDouble rho  = 0.001;
    //        fold(cellId, dist) = f(cellId, dist) + LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho); //slightly low
    //        if(m_normal[0] > 0) {
    //          // outlet
    //          if(dist == 5 || dist == 1 || dist == 5){
    //            fold(cellId, dist) = fold(cellId, LBMethod<LBTYPE>::oppositeDist(dist));
    //          }
    //          //          fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (1.0 - rho);
    //          //          if(dist == 4 || dist == 1 || dist == 5 || dist == 8) {
    //          //            fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
    //          //          } else
    //          //              if(dist == 7 || dist == 0 || dist == 6) {
    //          //            fold(cellId, dist) -= LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
    //          //          }
    //          //          if(dist == 8){
    //          fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (1.0 - rho);
    //          }
    //        } else {
    //          if(cellId == 341 || cellId == 0){
    //            fold(cellId, 8) += LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
    //            continue;
    //          }
    //          // inlet
    //          if(dist == 4 || dist == 1 || dist == 5 || dist == 8) {
    //            //          if(dist == 8) {
    ////            if((dist == 4) && cellId == 341){
    ////              continue;
    ////            }
    //            fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
    //            //            fold(cellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
    //            //          }
    //          } else if(dist == 0) {
    //            // outside direction
    //            //            fold(cellId, dist) += LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
    //            fold(cellId, dist) =
    //                -fold(cellId, LBMethod<LBTYPE>::oppositeDist(dist)) + LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
    //          } else if(dist == 7 || dist == 6 || dist == 3 || dist == 2) {
    ////            if((dist == 3 || dist==7) && cellId == 341){
    ////              continue;
    ////            }
    //            fold(cellId, dist) = LBMethod<LBTYPE>::m_weights[dist] * (targetRho - rho);
    //          }
    // if(dist == 7 || dist == 0 || dist == 6){
    //  fold(cellId, dist) = -fold(cellId, LBMethod<LBTYPE>::oppositeDist(dist)) + 0.001;
    //}
    //        }
    //      }
    //    }
    //    TERMM(-1, "TEsts");
  }


 private:
  GDouble              m_pressure = 1.0;
  VectorD<dim(LBTYPE)> m_normal;
  std::vector<GInt>    m_bndCells;
};

#endif // LBM_LBM_BND_IN_OUT_H
