#include "BCEuler2PhaseGhostPressure.h"
#include "EOS1Phase.h"

BCEuler2PhaseGhostPressure::BCEuler2PhaseGhostPressure(
  bool is_left,
  const EOS1Phase & eos_liq,
  const EOS1Phase & eos_vap,
  const FluxEuler2Phase & flux,
  const double & p)
  : BCEuler2PhaseGhost(is_left, eos_liq, eos_vap, flux),
    _p(p)
{
}

std::vector<double> BCEuler2PhaseGhostPressure::computeGhostSolution(
  const std::vector<double> & U, const double & A) const
{
  std::vector<double> U_ghost(7);

  const double aA_liq = U[0];
  const double a_liq = aA_liq / A;
  const double a_vap = 1.0 - a_liq;
  const double aA_vap = a_vap * A;

  const double r_liq = U[1] / aA_liq;
  const double u_liq = U[2] / U[1];
  const double e_liq = _eos_liq.e_from_p_r(_p, r_liq);
  const double E_liq = e_liq + 0.5 * u_liq * u_liq;

  const double r_vap = U[4] / aA_vap;
  const double u_vap = U[5] / U[4];
  const double e_vap = _eos_vap.e_from_p_r(_p, r_vap);
  const double E_vap = e_vap + 0.5 * u_vap * u_vap;

  U_ghost[0] = U[0];
  U_ghost[1] = U[1];
  U_ghost[2] = U[2];
  U_ghost[3] = aA_liq * r_liq * E_liq;
  U_ghost[4] = U[4];
  U_ghost[5] = U[5];
  U_ghost[6] = aA_vap * r_vap * E_vap;

  return U_ghost;
}
