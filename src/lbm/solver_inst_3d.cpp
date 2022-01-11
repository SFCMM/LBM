#include "solver.cpp"

// Main directions + diagonals (unsuitable for high Reynolds numbers)
extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes>;
template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes>;
extern template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes>;
template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes>;
extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes>;
template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes>;

// extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q19, LBEquation::Poisson>;
// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q19, LBEquation::Poisson>;
// extern template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q19, LBEquation::Poisson>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q19, LBEquation::Poisson>;
// extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q19, LBEquation::Poisson>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q19, LBEquation::Poisson>;

// extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes_Poisson>;
// extern template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes_Poisson>;
// extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q19, LBEquation::Navier_Stokes_Poisson>;

// Main directions + diagonals + tridiagonals
extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes>;
template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes>;
extern template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes>;
template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes>;
extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes>;
template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes>;

// extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q27, LBEquation::Poisson>;
// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q27, LBEquation::Poisson>;
// extern template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q27, LBEquation::Poisson>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q27, LBEquation::Poisson>;
// extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q27, LBEquation::Poisson>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q27, LBEquation::Poisson>;
//
// extern template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::no_debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes_Poisson>;
// extern template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes_Poisson>;
// extern template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes_Poisson>;
// template class LBMSolver<Debug_Level::max_debug, LBMethodType::D3Q27, LBEquation::Navier_Stokes_Poisson>;