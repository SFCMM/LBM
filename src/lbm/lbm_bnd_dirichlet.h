#ifndef LBM_LBM_BND_DIRICHLET_H
#define LBM_LBM_BND_DIRICHLET_H
#include <json.h>
#include <sfcmm_common.h>
#include "analytical_solutions.h"
#include "common/surface.h"
#include "lbm_bnd_interface.h"
#include "lbm_constants.h"
#include "lbm_variables.h"

using json = nlohmann::json;

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquation EQ>
class LBMBnd_DirichletNEEM : public LBMBndInterface {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<EQ, NDIM>;

 public:
  // todo: allow setting specified variables
  LBMBnd_DirichletNEEM(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties)
    : m_value(config::opt_config_value(properties, "value", NAN)), m_normal(surf->normal(0)) {
    for(const GInt cellId : surf->getCellList()) {
      m_bndCells.emplace_back(cellId);
    }

    for(const auto bndCellId : m_bndCells) {
      GInt extrapolationDir = -1;
      for(GInt dist = 0; dist < cartesian::maxNoNghbrs<NDIM>(); ++dist) {
        if(surf->neighbor(bndCellId, dist) == INVALID_CELLID) {
          if(extrapolationDir < 0) {
            extrapolationDir = dist;
          } else {
            // we are at a corner so we have to use a diagonal neighbor
            if(extrapolationDir == 0 && dist == 2) {
              extrapolationDir = 6;
            }
            if(extrapolationDir == 0 && dist == 3) {
              extrapolationDir = 7;
            }
            if(extrapolationDir == 1 && dist == 3) {
              extrapolationDir = 4;
            }
            if(extrapolationDir == 1 && dist == 2) {
              extrapolationDir = 5;
            }
          }
        }
      }
      const GInt extrapolationCellId = surf->neighbor(bndCellId, cartesian::oppositeDir<NDIM>(extrapolationDir));
      m_extrapolationCellId.emplace_back(extrapolationCellId);
      if(extrapolationCellId == INVALID_CELLID) {
        cerr0 << "bndCellId " << bndCellId << " extrapolationDir " << extrapolationDir << std::endl;
        TERMM(-1, "No valid extrapolation cellId");
      }
    }

    const GBool setAnalyticalValue = config::opt_config_value(properties, "setAnalyticalValue", false);

    if(setAnalyticalValue) {
      if constexpr(NDIM == 1) {
        m_value = analytical::poisson::poissonCHAI08_1(surf->center(m_bndCells[0]))[0];
      }
    }

    if(EQ != LBEquation::Poisson) {
      TERMM(-1, "FIX ME");
    }
  }
  ~LBMBnd_DirichletNEEM() override = default;

  // deleted constructors not needed
  LBMBnd_DirichletNEEM(const LBMBnd_DirichletNEEM&) = delete;
  LBMBnd_DirichletNEEM(LBMBnd_DirichletNEEM&&)      = delete;
  auto operator=(const LBMBnd_DirichletNEEM&) -> LBMBnd_DirichletNEEM& = delete;
  auto operator=(LBMBnd_DirichletNEEM&&) -> LBMBnd_DirichletNEEM& = delete;

  void initCnd(const std::function<GDouble&(GInt, GInt)>& vars) override {
    for(const auto cellId : m_bndCells) {
      vars(cellId, VAR::electricPotential()) = m_value;
    }
  }


  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& vars) override {
    for(GInt id = 0; id < m_bndCells.size(); ++id) {
      const GInt cellId              = m_bndCells[id];
      const GInt extraPolationCellId = m_extrapolationCellId[id];

      for(GInt dist = 0; dist < NDIST - 1; ++dist) {
        const GDouble weight = LBMethod<LBTYPE>::m_weights[dist];
        fold(cellId, dist) =
            weight * m_value + fold(extraPolationCellId, dist) - weight * vars(extraPolationCellId, VAR::electricPotential());
      }
      fold(cellId, NDIST - 1) = (LBMethod<LBTYPE>::m_weights[NDIST - 1] - 1.0) * m_value + fold(extraPolationCellId, NDIST - 1)
                                - (LBMethod<LBTYPE>::m_weights[NDIST - 1] - 1.0) * vars(extraPolationCellId, VAR::electricPotential());
    }
  }


 private:
  GDouble           m_value = 1.0;
  VectorD<NDIM>     m_normal;
  std::vector<GInt> m_bndCells;
  std::vector<GInt> m_extrapolationCellId;
};
#endif // LBM_LBM_BND_DIRICHLET_H
