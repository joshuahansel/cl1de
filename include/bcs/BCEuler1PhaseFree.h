#ifndef BCEuler1PhaseFree_H
#define BCEuler1PhaseFree_H

#include "BCEuler1Phase.h"

class EOS1Phase;

class BCEuler1PhaseFree : public BCEuler1Phase
{
public:
  BCEuler1PhaseFree(bool is_left, const EOS1Phase & eos);

  virtual std::vector<double> computeFlux(
    const std::vector<double> & U, const double & A) const override;

protected:
  const EOS1Phase & _eos;
};

#endif /* BCEuler1PhaseFree_H */
