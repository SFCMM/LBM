#ifndef LBM_LBM_BND_DIRICHLET_H
#define LBM_LBM_BND_DIRICHLET_H
#include <json.h>
#include <sfcmm_common.h>
#include "common/surface.h"
#include "lbm_bnd_interface.h"
#include "lbm_constants.h"
#include "lbm_variables.h"

using json = nlohmann::json;

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE>
class LBMBnd_DirichletNEEM : public LBMBndInterface {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<LBEquation::Navier_Stokes, NDIM>;

 public:
  // todo: allow setting specified variables
  LBMBnd_DirichletNEEM(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties) : m_value(properties["value"]) {
    for(const GInt cellId : surf->getCellList()) {
      m_bndCells.emplace_back(cellId);
    }
    if(NDIM == 2) {
      TERMM(-1, "FIX ME");
    }
  }
  ~LBMBnd_DirichletNEEM() override = default;

  // deleted constructors not needed
  LBMBnd_DirichletNEEM(const LBMBnd_DirichletNEEM&) = delete;
  LBMBnd_DirichletNEEM(LBMBnd_DirichletNEEM&&)      = delete;
  auto operator=(const LBMBnd_DirichletNEEM&) -> LBMBnd_DirichletNEEM& = delete;
  auto operator=(LBMBnd_DirichletNEEM&&) -> LBMBnd_DirichletNEEM& = delete;

  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& vars) override {
    for(const auto cellId : m_bndCells) {
      GInt extraPolationCellId = 0;
      if(cellId == 0) {
        extraPolationCellId = 1;
      } else {
        extraPolationCellId = cellId - 1;
      }
      for(GInt dist = 0; dist < NDIST - 1; ++dist) {
        const GDouble weight = LBMethod<LBTYPE>::m_weights[dist];
        // todo: fix this with correct access to the variables
        fold(cellId, dist) = weight * m_value + fold(extraPolationCellId, dist) - weight * vars(extraPolationCellId, 0);
      }
      fold(cellId, NDIST - 1) = (LBMethod<LBTYPE>::m_weights[NDIST - 1] - 1.0) * m_value + fold(extraPolationCellId, NDIST - 1)
                                - (LBMethod<LBTYPE>::m_weights[NDIST - 1] - 1.0) * vars(extraPolationCellId, 0);
    }
  }


 private:
  GDouble           m_value = 1.0;
  std::vector<GInt> m_bndCells;
};
#endif // LBM_LBM_BND_DIRICHLET_H
