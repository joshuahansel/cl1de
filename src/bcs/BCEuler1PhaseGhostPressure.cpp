#include "BCEuler1PhaseGhostPressure.h"
#include "EOS1Phase.h"

BCEuler1PhaseGhostPressure::BCEuler1PhaseGhostPressure(
  bool is_left,
  const EOS1Phase & eos,
  const FluxEuler1Phase & flux,
  const double & p)
  : BCEuler1PhaseGhost(is_left, eos, flux),
    _p(p)
{
}

std::vector<double> BCEuler1PhaseGhostPressure::computeGhostSolution(
  const std::vector<double> & U, const double & A) const
{
  std::vector<double> U_ghost(3);

  // density and velocity come from the solution
  const double r = U[0] / A;
  const double u = U[1] / U[0];
  const double e = _eos.e_from_p_r(_p, r);
  const double E = e + 0.5 * u * u;

  U_ghost[0] = U[0];
  U_ghost[1] = U[1];
  U_ghost[2] = U[0] * E;

  return U_ghost;
}
