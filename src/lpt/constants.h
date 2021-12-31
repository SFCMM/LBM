#ifndef LPT_CONSTANTS_H
#define LPT_CONSTANTS_H

enum class LPTType { Normal, High };

enum class LPTInitCond { load_from_CSV, randomvol_pos, randomplane_pos };

static constexpr std::array<std::string_view, 3> LPTInitCondName = {"load_from_CSV", "randomvol_pos", "randomplane_pos"};

static constexpr auto LPTInitCond(const std::string_view condName) -> LPTInitCond {
  if(condName == "load_from_CSV") {
    return LPTInitCond::load_from_CSV;
  }
  if(condName == "randomvol_pos") {
    return LPTInitCond::randomvol_pos;
  }
  if(condName == "randomplane_pos") {
    return LPTInitCond::randomplane_pos;
  }
  TERMM(-1, "Invalid initial condition configuration!");
}
#endif // LPT_CONSTANTS_H
