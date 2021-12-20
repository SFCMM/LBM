#include "solver.cpp"

extern template class LPTSolver<Debug_Level::no_debug, 2, LPTType::Normal>;
template class LPTSolver<Debug_Level::no_debug, 2, LPTType::Normal>;
extern template class LPTSolver<Debug_Level::debug, 2, LPTType::Normal>;
template class LPTSolver<Debug_Level::debug, 2, LPTType::Normal>;
extern template class LPTSolver<Debug_Level::max_debug, 2, LPTType::Normal>;
template class LPTSolver<Debug_Level::max_debug, 2, LPTType::Normal>;