#ifndef LBM_BND_NEUMANN_H
#define LBM_BND_NEUMANN_H
#include <json.h>
#include <sfcmm_common.h>
#include "analytical_solutions.h"
#include "bnd_dirichlet.h"
#include "bnd_interface.h"
#include "common/surface.h"
#include "constants.h"
#include "variables.h"

using json = nlohmann::json;

// template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ>
// class LBMBnd_DirichletBB : public LBMBndInterface {
//  private:
//   static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
//   static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
//
//   static constexpr GInt NVAR = noVars<LBTYPE>(EQ);
//
//   using VAR = LBMVariables<EQ, NDIM>;
//
//  public:
//   LBMBnd_DirichletBB(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties)
//     : m_bnd(surf), m_value(config::required_config_value<NVAR>(properties, "value")) {}
//   ~LBMBnd_DirichletBB() override = default;
//
//   // deleted constructors not needed
//   LBMBnd_DirichletBB(const LBMBnd_DirichletBB&) = delete;
//   LBMBnd_DirichletBB(LBMBnd_DirichletBB&&)      = delete;
//   auto operator=(const LBMBnd_DirichletBB&) -> LBMBnd_DirichletBB& = delete;
//   auto operator=(LBMBnd_DirichletBB&&) -> LBMBnd_DirichletBB& = delete;
//
//   void initCnd(const std::function<GDouble&(GInt, GInt)>& vars) override {
//     for(const auto bndCellId : m_bnd->getCellList()) {
//       for(GInt dir = 0; dir < NDIM; ++dir) {
//         vars(bndCellId, VAR::velocity(dir)) = m_value[dir];
//       }
//     }
//   }
//
//   void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
//                 const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}
//
//   void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold,
//              const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {
//     for(const auto bndCellId : m_bnd->getCellList()) {
//       for(GInt dist = 0; dist < NDIST - 1; ++dist) {
//         // skip setting dists that are orthogonal to the boundary
//         if(m_bnd->neighbor(bndCellId, dist) == INVALID_CELLID
//            && inDirection<NDIM>(VectorD<NDIM>(m_bnd->normal_p(bndCellId)), LBMethod<LBTYPE>::m_dirs[dist])) {
//           // dist in reflected/opposite direction i.e. to the inside of the bnd
//           GInt oppositeDist = LBMethod<LBTYPE>::oppositeDist(dist);
//
//           // standard bounceback i.e. distribution that hits the wall is reflected to the opposite distribution direction
//           fold(bndCellId, oppositeDist) = fpre(bndCellId, dist);
//
//           // set dirichlet value
//           if(NVAR == 1 && EQ == LBEquationType::Poisson) {
//             // todo: calculate diffusivity
//             const GDouble            dt          = 0.00392158 / 8.0; // lvl11
//             static constexpr GDouble k           = 27.79;
//             const GDouble            diffusivity = -k * 0.335 * dt; // lvl11:2.0002 with k = 2 * 27.79 1.99
//             fold(bndCellId, oppositeDist) -= LBMethod<LBTYPE>::m_poissonWeights[dist] * diffusivity * m_value[0];
//             TERMM(-1, "this is incorrect!");
//           }
//           if(EQ == LBEquationType::Navier_Stokes) {
//             // todo: calculate density
//             const GDouble density = 1.0;
//             for(GInt dir = 0; dir < NDIM; ++dir) {
//               fold(bndCellId, oppositeDist) += 2.0 / lbm_cssq * LBMethod<LBTYPE>::m_weights[oppositeDist] * density
//                                                * LBMethod<LBTYPE>::m_dirs[oppositeDist][dir] * m_value[dir];
//             }
//           }
//         }
//       }
//     }
//   }
//
//  private:
//   const SurfaceInterface* m_bnd = nullptr;
//
//   VectorD<NVAR> m_value;
// };

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ>
class LBMBnd_DirichletNEEM;

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ>
class LBMBnd_NeumannNEEM : public LBMBnd_DirichletNEEM<DEBUG_LEVEL, LBTYPE, EQ> {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using VAR = LBMVariables<EQ, NDIM>;

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

  void initCnd(const std::function<GDouble&(GInt, GInt)>& vars) override {}


  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) override {
    if(EQ != LBEquationType::Poisson) {
      TERMM(-1, "Not implemented");
    }

    // update value to the calculated gradient value
    for(GUint id = 0; id < this->no_cells(); ++id) {
      const GInt extrapolationId = this->extrapolationCellId(id);

      const GInt ex2Id = this->neighbor(extrapolationId, this->extrapolationDir(id));
      this->value(id) = (4.0 * vars(extrapolationId, VAR::electricPotential()) - vars(ex2Id, VAR::electricPotential()) + m_gradValue) / 3.0;
    }
    LBMBnd_DirichletNEEM<DEBUG_LEVEL, LBTYPE, EQ>::apply(f, fold, feq, vars);
  }


 private:
  GDouble m_gradValue = 0.0;
};
#endif // LBM_BND_NEUMANN_H
