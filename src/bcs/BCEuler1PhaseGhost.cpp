#include "BCEuler1PhaseGhost.h"
#include "FluxEuler1Phase.h"

#include <iostream>

BCEuler1PhaseGhost::BCEuler1PhaseGhost(
  bool is_left,
  const EOS1Phase & eos,
  const FluxEuler1Phase & flux)
  : BCEuler1Phase(is_left),
    _eos(eos),
    _flux(flux)
{
}

std::vector<double> BCEuler1PhaseGhost::computeFlux(
  const std::vector<double> & U, const double & A) const
{
  const std::vector<double> U_ghost = computeGhostSolution(U, A);

  if (_is_left)
    return _flux.computeFlux(U_ghost, U, A);
  else
    return _flux.computeFlux(U, U_ghost, A);
}
