#include "ICsEuler2PhaseRUP.h"
#include "EOS1Phase.h"
#include "Function.h"

ICsEuler2PhaseRUP::ICsEuler2PhaseRUP(
  const Function & a_liq_fn,
  const Function & r_liq_fn,
  const Function & u_liq_fn,
  const Function & p_liq_fn,
  const Function & r_vap_fn,
  const Function & u_vap_fn,
  const Function & p_vap_fn,
  const EOS1Phase & eos_liq,
  const EOS1Phase & eos_vap)
  : ICsEuler2Phase(eos_liq, eos_vap),

    _a_liq_fn(a_liq_fn),
    _r_liq_fn(r_liq_fn),
    _u_liq_fn(u_liq_fn),
    _p_liq_fn(p_liq_fn),
    _r_vap_fn(r_vap_fn),
    _u_vap_fn(u_vap_fn),
    _p_vap_fn(p_vap_fn)
{
}

void ICsEuler2PhaseRUP::computeICs(
  const double & x,
  const double & A,
  double & aA_liq,
  double & arA_liq,
  double & aruA_liq,
  double & arEA_liq,
  double & arA_vap,
  double & aruA_vap,
  double & arEA_vap) const
{
  const double a_liq = _a_liq_fn.value(x, 0);
  const double r_liq = _r_liq_fn.value(x, 0);
  const double u_liq = _u_liq_fn.value(x, 0);
  const double p_liq = _p_liq_fn.value(x, 0);
  const double r_vap = _r_vap_fn.value(x, 0);
  const double u_vap = _u_vap_fn.value(x, 0);
  const double p_vap = _p_vap_fn.value(x, 0);

  const double a_vap = 1 - a_liq;
  const double e_liq = _eos_liq.e_from_p_r(p_liq, r_liq);
  const double E_liq = e_liq + 0.5 * u_liq * u_liq;
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
