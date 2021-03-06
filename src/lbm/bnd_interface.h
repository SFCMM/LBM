#ifndef LBM_BND_INTERFACE_H
#define LBM_BND_INTERFACE_H

#include "constants.h"

class LBMBndInterface {
 public:
  LBMBndInterface()          = default;
  virtual ~LBMBndInterface() = default;

  // deleted constructors not needed
  LBMBndInterface(const LBMBndInterface&) = delete;
  LBMBndInterface(LBMBndInterface&&)      = delete;
  auto operator=(const LBMBndInterface&) -> LBMBndInterface& = delete;
  auto operator=(LBMBndInterface&&) -> LBMBndInterface& = delete;

  virtual void initCnd(const std::function<GDouble&(GInt, GInt)>& vars)                                                  = 0;
  virtual void preApply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                        const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars) = 0;
  virtual void apply(const std::function<GDouble&(GInt, GInt)>& f, const std::function<GDouble&(GInt, GInt)>& fold,
                     const std::function<GDouble&(GInt, GInt)>& feq, const std::function<GDouble&(GInt, GInt)>& vars)    = 0;

 private:
};

// todo: remove?
template <LBMethodType LBTYPE>
class LBMBndCell {
 public:
  LBMBndCell(const GInt mappedCellId, const VectorD<dim(LBTYPE)>& normal) : m_normal(normal), m_mappedCellId(mappedCellId) {}
  virtual ~LBMBndCell() = default;

  virtual void init() {
    // calculate local bndry normal
  }

  [[nodiscard]] auto mapped() const -> GInt { return m_mappedCellId; }
  [[nodiscard]] auto normal() const -> const VectorD<dim(LBTYPE)>& { return m_normal; }

  // todo: fix me
  //  deleted constructors not needed
  //  LBMBndCell(const LBMBndCell&) = delete;
  //  LBMBndCell(LBMBndCell&&)      = delete;
  //  virtual auto operator=(const LBMBndCell&) -> LBMBndCell& = delete;
  //  virtual auto operator=(LBMBndCell&&) -> LBMBndCell& = delete;

 private:
  VectorD<dim(LBTYPE)> m_normal;
  GInt                 m_mappedCellId = -1;
};
#endif // LBM_BND_INTERFACE_H
