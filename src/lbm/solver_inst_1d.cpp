#include "solver.cpp"

// todo: no test case
// extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes>;
// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes>;
// extern template class LBMSolver<Debug_Level::debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes>;
// extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes>;

extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D1Q3, LBEquation::Poisson>;
template class LBMSolver<Debug_Level::no_debug, LBMethodType::D1Q3, LBEquation::Poisson>;
extern template class LBMSolver<Debug_Level::debug, LBMethodType::D1Q3, LBEquation::Poisson>;
template class LBMSolver<Debug_Level::debug, LBMethodType::D1Q3, LBEquation::Poisson>;
extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D1Q3, LBEquation::Poisson>;
template class LBMSolver<Debug_Level::max_debug, LBMethodType::D1Q3, LBEquation::Poisson>;

// todo: no test case
// extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes_Poisson>;
// extern template class LBMSolver<Debug_Level::debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes_Poisson>;
// extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D1Q3, LBEquation::Navier_Stokes_Poisson>;