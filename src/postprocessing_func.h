#ifndef LBM_POSTPROCESSING_FUNC_H
#define LBM_POSTPROCESSING_FUNC_H

namespace pp {
enum class FuncType { LINE };

inline auto getFuncType(const GString& typeId) -> FuncType {
  if(typeId == "line") {
    return FuncType::LINE;
  }
  TERMM(-1, "Invalid functype");
}
}; // namespace pp

template <GInt NDIM>
class PostprocessFunctionInterface {
 public:
  PostprocessFunctionInterface() = default;
  virtual ~PostprocessFunctionInterface() = default;

  virtual void init() = 0;
  virtual void execute() = 0;

};

#endif // LBM_POSTPROCESSING_FUNC_H
