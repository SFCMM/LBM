#include "solver.cpp"

// todo: probably not needed
// extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes>;
// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes>;
// extern template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes>;
// extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes>;

extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q5, LBEquationType::Poisson>;
template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q5, LBEquationType::Poisson>;
extern template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q5, LBEquationType::Poisson>;
template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q5, LBEquationType::Poisson>;
extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q5, LBEquationType::Poisson>;
template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q5, LBEquationType::Poisson>;

// todo: probably not needed
// extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes_Poisson>;
// extern template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes_Poisson>;
// extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q5, LBEquation::Navier_Stokes_Poisson>;

extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes>;
template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes>;
extern template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes>;
template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes>;
extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes>;
template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes>;

extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q9, LBEquationType::Poisson>;
template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q9, LBEquationType::Poisson>;
extern template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q9, LBEquationType::Poisson>;
template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q9, LBEquationType::Poisson>;
extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q9, LBEquationType::Poisson>;
template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q9, LBEquationType::Poisson>;

extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes_Poisson>;
template class LBMSolver<Debug_Level::no_debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes_Poisson>;
extern template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes_Poisson>;
template class LBMSolver<Debug_Level::debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes_Poisson>;
extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes_Poisson>;
template class LBMSolver<Debug_Level::max_debug, LBMethodType::D2Q9, LBEquationType::Navier_Stokes_Poisson>;