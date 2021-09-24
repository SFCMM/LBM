#ifndef GRIDGENERATOR_WEIGHTMETHOD_H
#define GRIDGENERATOR_WEIGHTMETHOD_H

enum class WeightMethodsTypes { uniform };

class WeightMethod {
 public:
  virtual auto weight(GInt id) -> GFloat = 0;

  WeightMethod()                    = default;
  virtual ~WeightMethod()           = default;
  WeightMethod(const WeightMethod&) = delete;
  WeightMethod(WeightMethod&&)      = delete;
  auto operator=(const WeightMethod&) -> WeightMethod& = delete;
  auto operator=(WeightMethod&&) -> WeightMethod& = delete;
};


class WeightUniform : public WeightMethod {
 public:
  WeightUniform()                     = default;
  ~WeightUniform() override           = default;
  WeightUniform(const WeightUniform&) = delete;
  WeightUniform(WeightUniform&&)      = delete;
  auto operator=(const WeightUniform&) -> WeightUniform& = delete;
  auto operator=(WeightUniform&&) -> WeightUniform& = delete;

  auto weight(GInt /*id*/) -> GFloat override { return 1.0; }
};

#endif // GRIDGENERATOR_WEIGHTMETHOD_H
