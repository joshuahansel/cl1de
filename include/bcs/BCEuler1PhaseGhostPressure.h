#ifndef BCEuler1PhaseGhostPressure_H
#define BCEuler1PhaseGhostPressure_H

#include "BCEuler1PhaseGhost.h"

class BCEuler1PhaseGhostPressure : public BCEuler1PhaseGhost
{
public:
  BCEuler1PhaseGhostPressure(
    bool is_left, const EOS1Phase & eos, const FluxEuler1Phase & flux, const double & p);

protected:
  virtual std::vector<double> computeGhostSolution(
    const std::vector<double> & U, const double & A) const override;

  const double & _p;
};

#endif /* BCEuler1PhaseGhostPressure_H */
