#ifndef ICsEuler2PhaseRUP_H
#define ICsEuler2PhaseRUP_H

#include "ICsEuler2Phase.h"

class Function;

class ICsEuler2PhaseRUP : public ICsEuler2Phase
{
public:
  ICsEuler2PhaseRUP(
    const Function & a_liq_fn,
    const Function & r_liq_fn,
    const Function & u_liq_fn,
    const Function & p_liq_fn,
    const Function & r_vap_fn,
    const Function & u_vap_fn,
    const Function & p_vap_fn,
    const EOS1Phase & eos_liq,
    const EOS1Phase & eos_vap);

  virtual void computeICs(
    const double & x,
    const double & A,
    double & aA_liq,
    double & arA_liq,
    double & aruA_liq,
    double & arEA_liq,
    double & arA_vap,
    double & aruA_vap,
    double & arEA_vap) const override;

protected:
  const Function & _a_liq_fn;
  const Function & _r_liq_fn;
  const Function & _u_liq_fn;
  const Function & _p_liq_fn;
  const Function & _r_vap_fn;
  const Function & _u_vap_fn;
  const Function & _p_vap_fn;
};

#endif /* ICsEuler2PhaseRUP_H */
