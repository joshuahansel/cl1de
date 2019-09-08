#include "BCEuler2PhaseGhostFree.h"
#include "EOS1Phase.h"

BCEuler2PhaseGhostFree::BCEuler2PhaseGhostFree(
  bool is_left,
  const EOS1Phase & eos_liq,
  const EOS1Phase & eos_vap,
  const FluxEuler2Phase & flux)
  : BCEuler2PhaseGhost(is_left, eos_liq, eos_vap, flux)
{
}

std::vector<double> BCEuler2PhaseGhostFree::computeGhostSolution(
  const std::vector<double> & U, const double & A) const
{
  std::vector<double> U_ghost(7);

  U_ghost[0] = U[0];
  U_ghost[1] = U[1];
  U_ghost[2] = U[2];
  U_ghost[3] = U[3];
  U_ghost[4] = U[4];
  U_ghost[5] = U[5];
  U_ghost[6] = U[6];

  return U_ghost;
}
