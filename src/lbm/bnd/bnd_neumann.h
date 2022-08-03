#ifndef LBM_BND_NEUMANN_H
#define LBM_BND_NEUMANN_H
#include <json.h>
#include "bnd_dirichlet.h"
#include "bnd_interface.h"
#include "common/surface.h"
#include "lbm/analytical_solutions.h"
#include "lbm/constants.h"
#include "lbm/variables.h"
#include "sfcmm_common.h"

using json = nlohmann::json;

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ>
class LBMBnd_NeumannNEEM : public LBMBnd_DirichletNEEM<DEBUG_LEVEL, LBTYPE, EQ> {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR                  = LBMVariables<EQ, NDIM>;
  static constexpr GInt NVAR = noVars<LBTYPE>(EQ);

 public:
  // todo: allow setting specified variables
  LBMBnd_NeumannNEEM(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties)
    : LBMBnd_DirichletNEEM<DEBUG_LEVEL, LBTYPE, EQ>(surf, properties) {}
  ~LBMBnd_NeumannNEEM() override = default;

  // deleted constructors not needed
  LBMBnd_NeumannNEEM(const LBMBnd_NeumannNEEM&) = delete;
  LBMBnd_NeumannNEEM(LBMBnd_NeumannNEEM&&)      = delete;
  auto operator=(const LBMBnd_NeumannNEEM&) -> LBMBnd_NeumannNEEM& = delete;
  auto operator=(LBMBnd_NeumannNEEM&&) -> LBMBnd_NeumannNEEM& = delete;

  void initCnd(const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}


  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) override {
    // todo: implement
    if(EQ != LBEquationType::Poisson) {
      TERMM(-1, "Not implemented");
    }

    // update value to the calculated gradient value
    for(GInt id = 0; id < this->no_cells(); ++id) {
      const GInt extrapolationId = this->extrapolationCellId(id);

      const GInt ex2Id = this->neighbor(extrapolationId, this->extrapolationDir(id));
      calcDensity<NDIM, NDIST, EQ>(ex2Id, fold, vars);

      VectorD<NVAR> temp;
      temp[0]         = (4.0 * vars(extrapolationId, VAR::electricPotential()) - vars(ex2Id, VAR::electricPotential()) + m_gradValue) / 3.0;
      this->value(id) = temp;
    }
    LBMBnd_DirichletNEEM<DEBUG_LEVEL, LBTYPE, EQ>::apply(f, fold, feq, vars);
  }


 private:
  GDouble m_gradValue = 0.0;
};
#endif // LBM_BND_NEUMANN_H
