#ifndef BCEuler2PhaseGhostPressure_H
#define BCEuler2PhaseGhostPressure_H

#include "BCEuler2PhaseGhost.h"

class BCEuler2PhaseGhostPressure : public BCEuler2PhaseGhost
{
public:
  BCEuler2PhaseGhostPressure(
    bool is_left, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap, const FluxEuler2Phase & flux,
    const double & p);

protected:
  virtual std::vector<double> computeGhostSolution(
    const std::vector<double> & U, const double & A) const override;

  const double _p;
};

#endif /* BCEuler2PhaseGhostPressure_H */
