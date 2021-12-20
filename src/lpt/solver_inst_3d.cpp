#include "solver.cpp"

extern template class LPTSolver<Debug_Level::no_debug, 3, LPTType::Normal>;
template class LPTSolver<Debug_Level::no_debug, 3, LPTType::Normal>;
extern template class LPTSolver<Debug_Level::debug, 3, LPTType::Normal>;
template class LPTSolver<Debug_Level::debug, 3, LPTType::Normal>;
extern template class LPTSolver<Debug_Level::max_debug, 3, LPTType::Normal>;
template class LPTSolver<Debug_Level::max_debug, 3, LPTType::Normal>;