#ifndef ICsEuler1PhaseRUP_H
#define ICsEuler1PhaseRUP_H

#include "ICsEuler1Phase.h"

class Function;

class ICsEuler1PhaseRUP : public ICsEuler1Phase
{
public:
  ICsEuler1PhaseRUP(
    const Function & r_fn, const Function & u_fn, const Function & p_fn, const EOS1Phase & eos);

  virtual void computeICs(const double & x, const double & A, double & rA, double & ruA, double & rEA) const override;

protected:
  const Function & _r_fn;
  const Function & _u_fn;
  const Function & _p_fn;
};

#endif /* ICsEuler1PhaseRUP_H */
