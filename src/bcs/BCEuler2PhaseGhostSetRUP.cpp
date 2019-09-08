#include "BCEuler2PhaseGhostSetRUP.h"
#include "EOS1Phase.h"

BCEuler2PhaseGhostSetRUP::BCEuler2PhaseGhostSetRUP(
  bool is_left,
  const EOS1Phase & eos_liq,
  const EOS1Phase & eos_vap,
  const FluxEuler2Phase & flux,
  const double & a_liq,
  const double & r_liq, const double & u_liq, const double & p_liq,
  const double & r_vap, const double & u_vap, const double & p_vap)
  : BCEuler2PhaseGhost(is_left, eos_liq, eos_vap, flux),
    _a_liq(a_liq),
    _r_liq(r_liq),
    _u_liq(u_liq),
    _p_liq(p_liq),
    _r_vap(r_vap),
    _u_vap(u_vap),
    _p_vap(p_vap)
{
}

std::vector<double> BCEuler2PhaseGhostSetRUP::computeGhostSolution(
  const std::vector<double> & /*U*/, const double & A) const
{
  std::vector<double> U_ghost(7);

  const double a_vap = 1.0 - _a_liq;
  const double aA_liq = _a_liq * A;
  const double arA_liq = aA_liq * _r_liq;
  const double aruA_liq = arA_liq * _u_liq;
  const double E_liq = _eos_liq.e_from_p_r(_p_liq, _r_liq) + 0.5 * _u_liq * _u_liq;
  const double arEA_liq = arA_liq * E_liq;
  const double arA_vap = a_vap * _r_vap * A;
  const double aruA_vap = arA_vap * _u_vap;
  const double E_vap = _eos_vap.e_from_p_r(_p_vap, _r_vap) + 0.5 * _u_vap * _u_vap;
  const double arEA_vap = arA_vap * E_vap;

  U_ghost[0] = aA_liq;
  U_ghost[1] = arA_liq;
  U_ghost[2] = aruA_liq;
  U_ghost[3] = arEA_liq;
  U_ghost[4] = arA_vap;
  U_ghost[5] = aruA_vap;
  U_ghost[6] = arEA_vap;

  return U_ghost;
}
