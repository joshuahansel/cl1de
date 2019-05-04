#ifndef ICsEuler1PhasePUT_H
#define ICsEuler1PhasePUT_H

#include "ICsEuler1Phase.h"

class Function;

class ICsEuler1PhasePUT : public ICsEuler1Phase
{
public:
  ICsEuler1PhasePUT(
    const Function & p_fn, const Function & u_fn, const Function & T_fn, const EOS1Phase & eos);

  virtual void computeICs(const double & x, const double & A, double & rA, double & ruA, double & rEA) const override;

protected:
  const Function & _p_fn;
  const Function & _u_fn;
  const Function & _T_fn;
};

#endif /* ICsEuler1PhasePUT_H */
