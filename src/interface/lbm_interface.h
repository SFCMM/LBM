#ifndef LBM_INTERFACE_H
#define LBM_INTERFACE_H

class LBMethodInterface {
 public:
  LBMethodInterface()          = default;
  virtual ~LBMethodInterface() = default;

  // deleted constructors not needed
  LBMethodInterface(const LBMethodInterface&) = delete;
  LBMethodInterface(LBMethodInterface&&)      = delete;
  auto operator=(const LBMethodInterface&) -> LBMethodInterface& = delete;
  auto operator=(LBMethodInterface&&) -> LBMethodInterface& = delete;

  virtual void collisionStep() = 0;
};
#endif // LBM_INTERFACE_H
