#include "BCEuler1PhaseGhostWall.h"
#include "EOS1Phase.h"

BCEuler1PhaseGhostWall::BCEuler1PhaseGhostWall(
  bool is_left,
  const EOS1Phase & eos,
  const FluxEuler1Phase & flux)
  : BCEuler1PhaseGhost(is_left, eos, flux)
{
}

std::vector<double> BCEuler1PhaseGhostWall::computeGhostSolution(
  const std::vector<double> & U, const double & A) const
{
  std::vector<double> U_ghost(3);

  U_ghost[0] = U[0];
  U_ghost[1] = -U[1];
  U_ghost[2] = U[2];

  return U_ghost;
}
