#include "ICsEuler1PhaseRUP.h"
#include "EOS1Phase.h"
#include "Function.h"

ICsEuler1PhaseRUP::ICsEuler1PhaseRUP(
  const Function & r_fn, const Function & u_fn, const Function & p_fn, const EOS1Phase & eos)
  : ICsEuler1Phase(eos),

    _r_fn(r_fn),
    _u_fn(u_fn),
    _p_fn(p_fn)
{
}

void ICsEuler1PhaseRUP::computeICs(const double & x, const double & A, double & rA, double & ruA, double & rEA) const
{
  const double r = _r_fn.value(x, 0);
  const double u = _u_fn.value(x, 0);
  const double p = _p_fn.value(x, 0);

  const double e = _eos.e_from_p_r(p, r);
  const double E = e + 0.5 * u * u;

  rA = r * A;
  ruA = rA * u;
  rEA = rA * E;
}
