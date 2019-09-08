#include "ReconstructorEuler2PhaseSlopePUT.h"
#include "EOS1Phase.h"

ReconstructorEuler2PhaseSlopePUT::ReconstructorEuler2PhaseSlopePUT(
  const DoFHandlerEuler2Phase & dof_handler, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap)
  : ReconstructorEuler2PhaseSlope(dof_handler, eos_liq, eos_vap)
{
}

void ReconstructorEuler2PhaseSlopePUT::computeSlopeVariables(
  const double & aA_liq,
  const double & arA_liq,
  const double & aruA_liq,
  const double & arEA_liq,
  const double & arA_vap,
  const double & aruA_vap,
  const double & arEA_vap,
  const double & A,
  double & a_liq,
  double & p_liq,
  double & u_liq,
  double & T_liq,
  double & p_vap,
  double & u_vap,
  double & T_vap) const
{
  a_liq = aA_liq / A;
  const double a_vap = 1.0 - a_liq;

  const double r_liq = arA_liq / aA_liq;
  u_liq = aruA_liq / arA_liq;
  const double E_liq = arEA_liq / arA_liq;
  const double e_liq = E_liq - 0.5 * u_liq * u_liq;
  p_liq = _eos_liq.p_from_r_e(r_liq, e_liq);
  T_liq = _eos_liq.T_from_r_e(r_liq, e_liq);

  const double r_vap = arA_vap / (a_vap * A);
  u_vap = aruA_vap / arA_vap;
  const double E_vap = arEA_vap / arA_vap;
  const double e_vap = E_vap - 0.5 * u_vap * u_vap;
  p_vap = _eos_vap.p_from_r_e(r_vap, e_vap);
  T_vap = _eos_vap.T_from_r_e(r_vap, e_vap);
}

void ReconstructorEuler2PhaseSlopePUT::computeConservativeVariables(
  const double & a_liq,
  const double & p_liq,
  const double & u_liq,
  const double & T_liq,
  const double & p_vap,
  const double & u_vap,
  const double & T_vap,
  const double & A,
  double & aA_liq,
  double & arA_liq,
  double & aruA_liq,
  double & arEA_liq,
  double & arA_vap,
  double & aruA_vap,
  double & arEA_vap) const
{
  const double a_vap = 1.0 - a_liq;
  const double r_liq = _eos_liq.r_from_p_T(p_liq, T_liq);
  const double e_liq = _eos_liq.e_from_p_r(p_liq, r_liq);
  const double E_liq = e_liq + 0.5 * u_liq * u_liq;
  const double r_vap = _eos_vap.r_from_p_T(p_vap, T_vap);
  const double e_vap = _eos_vap.e_from_p_r(p_vap, r_vap);
  const double E_vap = e_vap + 0.5 * u_vap * u_vap;

  aA_liq = a_liq * A;
  arA_liq = aA_liq * r_liq;
  aruA_liq = arA_liq * u_liq;
  arEA_liq = arA_liq * E_liq;
  arA_vap = a_vap * r_vap * A;
  aruA_vap = arA_vap * u_vap;
  arEA_vap = arA_vap * E_vap;
}
