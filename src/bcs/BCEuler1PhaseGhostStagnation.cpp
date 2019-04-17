#include "BCEuler1PhaseGhostStagnation.h"
#include "EOS1Phase.h"

BCEuler1PhaseGhostStagnation::BCEuler1PhaseGhostStagnation(
  bool is_left,
  const EOS1Phase & eos,
  const FluxEuler1Phase & flux,
  const double & p0,
  const double & T0)
  : BCEuler1PhaseGhost(is_left, eos, flux),
    _p0(p0),
    _T0(T0)
{
}

std::vector<double> BCEuler1PhaseGhostStagnation::computeGhostSolution(
  const std::vector<double> & U, const double & A) const
{
  std::vector<double> U_ghost(3);

  // velocity comes from the solution
  const double u = U[1] / U[0];
  const double h0 = _eos.h_from_p_T(_p0, _T0);
  const double h = h0 - 0.5 * u * u;
  const double s = _eos.s_from_p_T(_p0, _T0);
  const double r = _eos.r_from_h_s(h, s);
  const double e = _eos.e_from_r_h(r, h);
  const double E = e + 0.5 * u * u;

  U_ghost[0] = r * A;
  U_ghost[1] = U_ghost[0] * u;
  U_ghost[2] = U_ghost[0] * E;

  return U_ghost;
}
