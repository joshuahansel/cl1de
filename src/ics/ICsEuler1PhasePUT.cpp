#include "ICsEuler1PhasePUT.h"
#include "EOS1Phase.h"
#include "Function.h"

ICsEuler1PhasePUT::ICsEuler1PhasePUT(
  const Function & p_fn, const Function & u_fn, const Function & T_fn, const EOS1Phase & eos)
  : ICsEuler1Phase(eos),

    _p_fn(p_fn),
    _u_fn(u_fn),
    _T_fn(T_fn)
{
}

void ICsEuler1PhasePUT::computeICs(const double & x, const double & A, double & rA, double & ruA, double & rEA) const
{
  const double p = _p_fn.value(x, 0);
  const double u = _u_fn.value(x, 0);
  const double T = _T_fn.value(x, 0);

  const double r = _eos.r_from_p_T(p, T);
  const double e = _eos.e_from_p_r(p, r);
  const double E = e + 0.5 * u * u;

  rA = r * A;
  ruA = rA * u;
  rEA = rA * E;
}
