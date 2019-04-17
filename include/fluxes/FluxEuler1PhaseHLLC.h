#ifndef FluxEuler1PhaseHLLC_H
#define FluxEuler1PhaseHLLC_H

#include "FluxEuler1Phase.h"

class FluxEuler1PhaseHLLC : public FluxEuler1Phase
{
public:
  FluxEuler1PhaseHLLC(const EOS1Phase & eos);

  virtual std::vector<double> computeFlux(
    const std::vector<double> & U_L,
    const std::vector<double> & U_R,
    const double & A) const override;
};

#endif /* FluxEuler1PhaseHLLC_H */
