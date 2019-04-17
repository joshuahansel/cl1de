#ifndef FluxEuler1PhaseCentered_H
#define FluxEuler1PhaseCentered_H

#include "FluxEuler1Phase.h"

class FluxEuler1PhaseCentered : public FluxEuler1Phase
{
public:
  FluxEuler1PhaseCentered(const EOS1Phase & eos);

  virtual std::vector<double> computeFlux(
    const std::vector<double> & U_L,
    const std::vector<double> & U_R,
    const double & A) const override;
};

#endif /* FluxEuler1PhaseCentered_H */
