#ifndef BCEuler1PhaseGhostStagnation_H
#define BCEuler1PhaseGhostStagnation_H

#include "BCEuler1PhaseGhost.h"

class BCEuler1PhaseGhostStagnation : public BCEuler1PhaseGhost
{
public:
  BCEuler1PhaseGhostStagnation(
    bool is_left,
    const EOS1Phase & eos,
    const FluxEuler1Phase & flux,
    const double & p0,
    const double & T0);

protected:
  virtual std::vector<double> computeGhostSolution(
    const std::vector<double> & U, const double & A) const override;

  const double & _p0;
  const double & _T0;
};

#endif /* BCEuler1PhaseGhostStagnation_H */
