#ifndef LBM_MATHEXPR_H
#define LBM_MATHEXPR_H
#include <utility>
#include <sfcmm_common.h>
#include "exprtk.h"

template <GInt NDIM>
class MathExpression {
 public:
  MathExpression() = default;
  MathExpression(GString expr_string) : m_expr_str(std::move(expr_string)) {
    m_sym.add_variable("x", m_params[0]);
    if constexpr(NDIM > 1) {
      m_sym.add_variable("y", m_params[1]);
      if constexpr(NDIM == 3) {
        m_sym.add_variable("z", m_params[2]);
      }
    }
    m_sym.add_constants();

    m_expr.register_symbol_table(m_sym);
    m_parse.compile(m_expr_str, m_expr);
  }

  auto eval(const std::array<GDouble, NDIM>& x) -> GDouble {
    m_params = x;
    return m_expr.value();
  }

  auto eval(const VectorD<NDIM>& x) -> GDouble {
    for(GInt dir = 0; dir < NDIM; ++dir) {
      m_params[dir] = x[dir];
    }
    return m_expr.value();
  }

  auto empty() -> GBool { return m_expr_str.empty(); }

 private:
  GString                   m_expr_str;
  std::array<GDouble, NDIM> m_params;

  // storage for exprtk
  exprtk::symbol_table<GDouble> m_sym;
  exprtk::expression<GDouble>   m_expr;
  exprtk::parser<GDouble>       m_parse;
};

#endif // LBM_MATHEXPR_H
