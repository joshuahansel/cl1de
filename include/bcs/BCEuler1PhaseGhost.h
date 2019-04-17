#ifndef BCEuler1PhaseGhost_H
#define BCEuler1PhaseGhost_H

#include "BCEuler1Phase.h"

class EOS1Phase;
class FluxEuler1Phase;

class BCEuler1PhaseGhost : public BCEuler1Phase
{
public:
  BCEuler1PhaseGhost(bool is_left, const EOS1Phase & eos, const FluxEuler1Phase & flux);

  virtual std::vector<double> computeFlux(
    const std::vector<double> & U, const double & A) const override;

protected:
  virtual std::vector<double> computeGhostSolution(
    const std::vector<double> & U, const double & A) const = 0;

  const EOS1Phase & _eos;
  const FluxEuler1Phase & _flux;
};

#endif /* BCEuler1PhaseGhost_H */
