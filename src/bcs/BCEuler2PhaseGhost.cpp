#include "BCEuler2PhaseGhost.h"
#include "FluxEuler2Phase.h"
#include "UtilsEuler2Phase.h"

BCEuler2PhaseGhost::BCEuler2PhaseGhost(
  bool is_left,
  const EOS1Phase & eos_liq,
  const EOS1Phase & eos_vap,
  const FluxEuler2Phase & flux)
  : BCEuler2Phase(is_left),
    _eos_liq(eos_liq),
    _eos_vap(eos_vap),
    _flux(flux)
{
}

std::vector<double> BCEuler2PhaseGhost::computeFlux(
  const std::vector<double> & U, const double & A) const
{
  const std::vector<double> U_ghost = computeGhostSolution(U, A);

  std::vector<double> f_L(Euler2Phase::n_eq);
  std::vector<double> f_R(Euler2Phase::n_eq);
  if (_is_left)
  {
    _flux.computeFlux(U_ghost, U, A, f_L, f_R);
    return f_R;
  }
  else
  {
    _flux.computeFlux(U, U_ghost, A, f_L, f_R);
    return f_L;
  }
}
