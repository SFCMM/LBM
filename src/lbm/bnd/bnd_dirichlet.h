#ifndef LBM_BND_DIRICHLET_H
#define LBM_BND_DIRICHLET_H
#include <json.h>
#include "bnd_interface.h"
#include "bnd_wetnode.h"
#include "common/surface.h"
#include "lbm/analytical_solutions.h"
#include "lbm/constants.h"
#include "lbm/equilibrium_func.h"
#include "lbm/moments.h"
#include "lbm/variables.h"
#include "sfcmm_common.h"

using json = nlohmann::json;

/// Dirichlet Boundary condition implemented using the Bounce-back method
/// \tparam DEBUG_LEVEL
/// \tparam LBTYPE
/// \tparam EQ
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ>
class LBMBnd_DirichletBB : public LBMBndInterface {
 private:
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;

  static constexpr GInt NVAR = noVars<LBTYPE>(EQ);

  using VAR = LBMVariables<EQ, NDIM>;

 public:
  LBMBnd_DirichletBB(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties)
    : m_bnd(surf), m_value(config::required_config_value<NVAR>(properties, "value")) {}
  LBMBnd_DirichletBB(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& /*properties*/, const GDouble value) : m_bnd(surf) {
    m_value.fill(value);
  }
  ~LBMBnd_DirichletBB() override = default;

  // deleted constructors not needed
  LBMBnd_DirichletBB(const LBMBnd_DirichletBB&)                    = delete;
  LBMBnd_DirichletBB(LBMBnd_DirichletBB&&)                         = delete;
  auto operator=(const LBMBnd_DirichletBB&) -> LBMBnd_DirichletBB& = delete;
  auto operator=(LBMBnd_DirichletBB&&) -> LBMBnd_DirichletBB&      = delete;

  void initCnd(const std::function<GDouble&(GInt, GInt)>& vars) override {
    for(const auto bndCellId : m_bnd->getCellList()) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(bndCellId, VAR::velocity(dir)) = m_value[dir];
      }
    }
  }

  /// Nothing to do here
  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  /// Apply the Dirichlet boundary condition (setting a constant value for the result)
  /// \param fpre
  /// \param fold
  /// \param feq
  /// \param vars
  void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) override {
    // todo: this is not correct
    for(const auto bndCellId : m_bnd->getCellList()) {
      // using the general formulation for bb dirichlet bndcnd
      apply(bndCellId, m_value.data(), fpre, fold, feq, vars);
    }
  }


 protected:
  /// Implementation of the Dirichlet boundary condition.
  /// \tparam VALZERO Set 0 which means setting only bounce-back for this bnd
  /// \tparam SCALAR We are setting a scalar value
  /// \param bndCellId Id of the cell to be set
  /// \param values Dirichlet value to be set
  /// \param fpre Distribution
  /// \param fold Distribution
  template <GBool VALZERO = false, GBool SCALAR = false>
  void apply(const GInt bndCellId, const GDouble* values, const std::function<GDouble&(GInt, GInt)>& fpre,
             const std::function<GDouble&(GInt, GInt)>& fold, const std::function<GDouble&(GInt, GInt)>& /*feq*/,
             const std::function<GDouble&(GInt, GInt)>& /*vars*/) {
    for(GInt dist = 0; dist < NDIST - 1; ++dist) {
      // skip setting dists that are orthogonal to the boundary
      if(m_bnd->neighbor(bndCellId, dist) == INVALID_CELLID
         && inDirection<NDIM>(VectorD<NDIM>(m_bnd->normal_p(bndCellId)), LBMethod<LBTYPE>::m_dirs[dist])) {
        // dist in reflected/opposite direction i.e. to the inside of the bnd
        GInt oppositeDist = LBMethod<LBTYPE>::oppositeDist(dist);

        // standard bounceback i.e. distribution that hits the wall is reflected to the opposite distribution direction
        fold(bndCellId, oppositeDist) = fpre(bndCellId, dist);

        // todo: move specialization to a extra function

        // set dirichlet value for the poisson equation
        if(NVAR == 1 && EQ == LBEquationType::Poisson && !VALZERO) {
          // todo: calculate diffusivity
          const GDouble            dt          = 0.00392158 / 8.0; // lvl11
          static constexpr GDouble k           = 27.79;
          const GDouble            diffusivity = -k * 0.335 * dt; // lvl11:2.0002 with k = 2 * 27.79 1.99
          fold(bndCellId, oppositeDist) -= LBMethod<LBTYPE>::m_poissonWeights[dist] * diffusivity * m_value[0];
          TERMM(-1, "this is incorrect!");
        }

        // set dirichlet value for the Navier stokes equation
        if(EQ == LBEquationType::Navier_Stokes && !VALZERO) {
          // todo: calculate density
          // todo: cleanup
          const GDouble density = 1.0;
          if(SCALAR) {
            fold(bndCellId, oppositeDist) += density * 2.0 / lbm_cssq * LBMethod<LBTYPE>::m_weights[oppositeDist] * values[oppositeDist];
          } else {
            for(GInt dir = 0; dir < NDIM; ++dir) {
              fold(bndCellId, oppositeDist) += density * 2.0 / lbm_cssq * LBMethod<LBTYPE>::m_weights[oppositeDist]
                                               * LBMethod<LBTYPE>::m_dirs[oppositeDist][dir] * values[dir];
            }
          }
        }
      }
    }
  }

 private:
  const SurfaceInterface* m_bnd = nullptr;

  VectorD<NVAR> m_value;
};

/// Dirichlet Boundary condition implemented using the Equilibrium method
/// \tparam DEBUG_LEVEL
/// \tparam LBTYPE
/// \tparam EQ
template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ>
class LBMBnd_DirichletEQ : public LBMBndInterface, protected LBMBnd_wallWetnode<DEBUG_LEVEL, LBTYPE> {
 private:
  using method = LBMethod<LBTYPE>;

  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;

  static constexpr GInt NVAR = noVars<LBTYPE>(EQ);

  using VAR = LBMVariables<EQ, NDIM>;

 public:
  LBMBnd_DirichletEQ(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties)
    : LBMBnd_wallWetnode<DEBUG_LEVEL, LBTYPE>(surf), m_bnd(surf), m_value(config::required_config_value<NVAR>(properties, "value")) {
    init();
  }
  LBMBnd_DirichletEQ(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& /*properties*/, const GDouble value)
    : LBMBnd_wallWetnode<DEBUG_LEVEL, LBTYPE>(surf), m_bnd(surf) {
    m_value.fill(value);
    init();
  }
  LBMBnd_DirichletEQ(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const GDouble value)
    : LBMBnd_wallWetnode<DEBUG_LEVEL, LBTYPE>(surf), m_bnd(surf) {
    m_value.fill(value);
    init();
  }
  ~LBMBnd_DirichletEQ() override = default;

  // deleted constructors not needed
  LBMBnd_DirichletEQ(const LBMBnd_DirichletEQ&)                    = delete;
  LBMBnd_DirichletEQ(LBMBnd_DirichletEQ&&)                         = delete;
  auto operator=(const LBMBnd_DirichletEQ&) -> LBMBnd_DirichletEQ& = delete;
  auto operator=(LBMBnd_DirichletEQ&&) -> LBMBnd_DirichletEQ&      = delete;

  virtual void init() {
    GInt index = 0;
    for(const auto bndCellId : m_bnd->getCellList()) {
      cell2Bnd[bndCellId] = index;
      ++index;
    }
  }

  void initCnd(const std::function<GDouble&(GInt, GInt)>& vars) override {
    // todo: this is not correct since it only allows setting the velocity
    for(const auto bndCellId : m_bnd->getCellList()) {
      for(GInt dir = 0; dir < NDIM; ++dir) {
        vars(bndCellId, VAR::velocity(dir)) = m_value[dir];
      }
    }
  }

  /// Nothing to do here
  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  /// Apply the Dirichlet boundary condition (setting a constant value for the result)
  /// \param fpre
  /// \param fold
  /// \param feq
  /// \param vars
  void apply(const std::function<GDouble&(GInt, GInt)>& fpre, const std::function<GDouble&(GInt, GInt)>& fold,
             const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) override {
    for(const auto bndCellId : m_bnd->getCellList()) {
      // using the general formulation for bb dirichlet bndcnd
      apply(bndCellId, m_value.data(), fpre, fold, feq, vars);
    }
  }


 protected:
  /// Implementation of the Dirichlet boundary condition.
  /// \tparam VALZERO Set 0 which means setting only the qauilibrium value for this bnd
  /// \tparam SCALAR We are setting a scalar value
  /// \param bndCellId Id of the cell to be set
  /// \param values Dirichlet value to be set
  /// \param fpre Distribution
  /// \param fold Distribution
  template <GBool VALZERO = false>
  void apply(const GInt bndCellId, const GDouble* values, const std::function<GDouble&(GInt, GInt)>& /*fpre*/,
             const std::function<GDouble&(GInt, GInt)>& fold, const std::function<GDouble&(GInt, GInt)>& /*feq*/,
             const std::function<GDouble&(GInt, GInt)>& vars) {
    const GInt index = cell2Bnd[bndCellId];
    // at wall set wall velocity
    // todo: allow setting the variable to be set
    for(GInt dir = 0; dir < NDIM; ++dir) {
      vars(bndCellId, VAR::velocity(dir)) = values[dir];
    }


    // update the density value with velocity set to 0
    if(VALZERO) {
      calcDensity_limited<NDIM, NDIST, LBEquationType::Navier_Stokes, true>(bndCellId, this->limitedDist()[index],
                                                                            this->limitedConst()[index], nullptr, fold, vars);
    } else {
      calcDensity_limited<NDIM, NDIST, LBEquationType::Navier_Stokes, false>(
          bndCellId, this->limitedDist()[index], this->limitedConst()[index], &this->normal()[index][0], fold, vars);
    }


    for(GInt dist = 0; dist < NDIST; ++dist) {
      // todo: make settable
      if(VALZERO) {
        fold(bndCellId, dist) = eq::defaultEq(method::m_weights[dist], vars(bndCellId, VAR::rho()));
      } else {
        eq::defaultEq<LBTYPE>(&fold(bndCellId, 0), vars(bndCellId, VAR::rho()), &vars(bndCellId, VAR::velocity(0)));
      }
    }
  }

 private:
  const SurfaceInterface* m_bnd = nullptr;

  VectorD<NVAR>                  m_value;
  std::unordered_map<GInt, GInt> cell2Bnd;
};

template <Debug_Level DEBUG_LEVEL, LBMethodType LBTYPE, LBEquationType EQ>
class LBMBnd_DirichletNEEM : public LBMBndInterface {
 private:
  using method                = LBMethod<LBTYPE>;
  static constexpr GInt NDIM  = LBMethod<LBTYPE>::m_dim;
  static constexpr GInt NDIST = LBMethod<LBTYPE>::m_noDists;

  using METH                 = LBMethod<getLBMethodType(NDIM, NDIST)>;
  using VAR                  = LBMVariables<EQ, NDIM>;
  static constexpr GInt NVAR = noVars<LBTYPE>(EQ);

 public:
  // todo: allow setting specified variables
  LBMBnd_DirichletNEEM(const Surface<DEBUG_LEVEL, dim(LBTYPE)>* surf, const json& properties) : m_normal(surf->normal()), m_bnd(surf) {
    m_value.resize(surf->no_cells());

    // todo: move this to a function and object
    if(properties["value"].is_string() || (properties["value"].is_array() && properties["value"][0].is_string())) {
      if(properties["value"].is_array()) {
        TERMM(-1, "Impl");
      } else {
        m_mathExpression[0] = config::required_math_expression<NDIM>(properties, "value");
      }
      GInt bndId = 0;
      for(const auto cellId : m_bnd->getCellList()) {
        for(GInt dir = 0; dir < NVAR; ++dir) {
          m_value[bndId][dir] = m_mathExpression[dir]->eval(surf->center(cellId));
        }
        ++bndId;
      }
    } else {
      const VectorD<NVAR> value_to_set = (config::required_config_value<NVAR>(properties, "value"));
      std::fill(m_value.begin(), m_value.end(), value_to_set);
    }

    for(const auto bndCellId : m_bnd->getCellList()) {
      GInt extrapolationDir = -1;
      for(GInt dist = 0; dist < cartesian::maxNoNghbrs<NDIM>(); ++dist) {
        if(surf->neighbor(bndCellId, dist) == INVALID_CELLID) {
          if(extrapolationDir < 0) {
            extrapolationDir = dist;
          } else {
            // todo: move to function
            //  we are at a corner so we have to use a diagonal neighbor
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
      m_extrapolationDir.emplace_back(cartesian::oppositeDir<NDIM>(extrapolationDir));
      if(extrapolationCellId == INVALID_CELLID) {
        cerr0 << "bndCellId " << bndCellId << " extrapolationDir " << extrapolationDir << std::endl;
        TERMM(-1, "No valid extrapolation cellId");
      }
    }

    const GBool setAnalyticalValue = config::opt_config_value(properties, "setAnalyticalValue", false);

    // todo: move this somewhere else
    if(setAnalyticalValue) {
      if constexpr(NDIM == 1) {
        GInt bndId = 0;
        for(const auto bndCellId : m_bnd->getCellList()) {
          m_value[bndId++] = analytical::poisson::poissonCHAI08_1(surf->center(bndCellId));
        }
      }
    }

    if(EQ != LBEquationType::Poisson) {
      TERMM(-1, "FIX ME");
    }
  }
  ~LBMBnd_DirichletNEEM() override = default;

  // deleted constructors not needed
  LBMBnd_DirichletNEEM(const LBMBnd_DirichletNEEM&)                    = delete;
  LBMBnd_DirichletNEEM(LBMBnd_DirichletNEEM&&)                         = delete;
  auto operator=(const LBMBnd_DirichletNEEM&) -> LBMBnd_DirichletNEEM& = delete;
  auto operator=(LBMBnd_DirichletNEEM&&) -> LBMBnd_DirichletNEEM&      = delete;

  void initCnd(const std::function<GDouble&(GInt, GInt)>& vars) override {
    GInt bndId = 0;
    for(const auto cellId : m_bnd->getCellList()) {
      vars(cellId, VAR::electricPotential()) = m_value[bndId][0];
      ++bndId;
    }
  }


  void preApply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>& /*fold*/,
                const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& /*vars*/) override {}

  void apply(const std::function<GDouble&(GInt, GInt)>& /*f*/, const std::function<GDouble&(GInt, GInt)>&   fold,
             const std::function<GDouble&(GInt, GInt)>& /*feq*/, const std::function<GDouble&(GInt, GInt)>& vars) override {
    // recalculate density because fold is post streaming!
    calcDensity<NDIM, NDIST, EQ>(m_extrapolationCellId, fold, vars);

    GInt bndId = 0;
    for(const auto cellId : m_bnd->getCellList()) {
      const GInt extraPolationCellId         = m_extrapolationCellId[bndId];
      vars(cellId, VAR::electricPotential()) = m_value[bndId][0];

      for(GInt dist = 0; dist < NDIST - 1; ++dist) {
        const GDouble weight = LBMethod<LBTYPE>::m_weights[dist];
        fold(cellId, dist) =
            weight * m_value[bndId][0] + fold(extraPolationCellId, dist) - weight * vars(extraPolationCellId, VAR::electricPotential());
      }
      fold(cellId, NDIST - 1) = (LBMethod<LBTYPE>::m_weights[NDIST - 1] - 1.0) * m_value[bndId][0] + fold(extraPolationCellId, NDIST - 1)
                                - (LBMethod<LBTYPE>::m_weights[NDIST - 1] - 1.0) * vars(extraPolationCellId, VAR::electricPotential());
      ++bndId;
    }
  }

 protected:
  [[nodiscard]] auto no_cells() const -> GInt { return m_bnd->no_cells(); }

  [[nodiscard]] auto extrapolationCellId(const GInt id) const -> GInt { return m_extrapolationCellId[id]; }

  [[nodiscard]] auto extrapolationDir(const GInt id) const -> GInt { return m_extrapolationDir[id]; }

  [[nodiscard]] auto neighbor(const GInt cellId, const GInt dir) const -> GInt { return m_bnd->neighbor(cellId, dir); }

  auto value(const GInt id) -> VectorD<NVAR>& { return m_value[id]; }


 private:
  std::vector<VectorD<NVAR>>                              m_value;
  std::array<std::unique_ptr<MathExpression<NDIM>>, NDIM> m_mathExpression;
  VectorD<NDIM>                                           m_normal;
  std::vector<GInt>                                       m_extrapolationCellId;
  std::vector<GInt>                                       m_extrapolationDir;
  const SurfaceInterface*                                 m_bnd = nullptr;
};
#endif // LBM_BND_DIRICHLET_H
