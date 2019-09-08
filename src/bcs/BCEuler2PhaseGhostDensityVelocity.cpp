#include "BCEuler2PhaseGhostDensityVelocity.h"
#include "EOS1Phase.h"

BCEuler2PhaseGhostDensityVelocity::BCEuler2PhaseGhostDensityVelocity(
  bool is_left,
  const EOS1Phase & eos_liq,
  const EOS1Phase & eos_vap,
  const FluxEuler2Phase & flux,
  const double & a_liq,
  const double & r_liq, const double & u_liq,
  const double & r_vap, const double & u_vap)
  : BCEuler2PhaseGhost(is_left, eos_liq, eos_vap, flux),
    _a_liq(a_liq),
    _r_liq(r_liq),
    _u_liq(u_liq),
    _r_vap(r_vap),
    _u_vap(u_vap)
{
}

std::vector<double> BCEuler2PhaseGhostDensityVelocity::computeGhostSolution(
  const std::vector<double> & U, const double & A) const
{
  std::vector<double> U_ghost(7);

  const double aA_liq   = U[0];
  const double arA_liq  = U[1];
  const double aruA_liq = U[2];
  const double arEA_liq = U[3];
  const double arA_vap  = U[4];
  const double aruA_vap = U[5];
  const double arEA_vap = U[6];

  const double a_liq = aA_liq / A;
  const double a_vap = 1.0 - a_liq;
  const double aA_vap = a_vap * A;
  const double aA_liq_b = _a_liq * A;
  const double aA_vap_b = (1.0 - _a_liq) * A;

  const double r_liq = arA_liq / aA_liq;
  const double u_liq = aruA_liq / arA_liq;
  const double E_liq = arEA_liq / arA_liq;
  const double e_liq = E_liq - 0.5 * u_liq * u_liq;
  const double p_liq = _eos_liq.p_from_r_e(r_liq, e_liq);
  const double e_liq_b = _eos_liq.e_from_p_r(p_liq, _r_liq);
  const double E_liq_b = e_liq_b + 0.5 * _u_liq * _u_liq;
  const double arA_liq_b = aA_liq_b * _r_liq;
  const double aruA_liq_b = arA_liq_b * _u_liq;
  const double arEA_liq_b = arA_liq_b * E_liq_b;

  const double r_vap = arA_vap / aA_vap;
  const double u_vap = aruA_vap / arA_vap;
  const double E_vap = arEA_vap / arA_vap;
  const double e_vap = E_vap - 0.5 * u_vap * u_vap;
  const double p_vap = _eos_vap.p_from_r_e(r_vap, e_vap);
  const double e_vap_b = _eos_vap.e_from_p_r(p_vap, _r_vap);
  const double E_vap_b = e_vap_b + 0.5 * _u_vap * _u_vap;
  const double arA_vap_b = aA_vap_b * _r_vap;
  const double aruA_vap_b = arA_vap_b * _u_vap;
  const double arEA_vap_b = arA_vap_b * E_vap_b;

  U_ghost[0] = aA_liq_b;
  U_ghost[1] = arA_liq_b;
  U_ghost[2] = aruA_liq_b;
  U_ghost[3] = arEA_liq_b;
  U_ghost[4] = arA_vap_b;
  U_ghost[5] = aruA_vap_b;
  U_ghost[6] = arEA_vap_b;

  return U_ghost;
}
