#ifndef BCEuler2PhaseGhostFree_H
#define BCEuler2PhaseGhostFree_H

#include "BCEuler2PhaseGhost.h"

class BCEuler2PhaseGhostFree : public BCEuler2PhaseGhost
{
public:
  BCEuler2PhaseGhostFree(
    bool is_left, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap, const FluxEuler2Phase & flux);

protected:
  virtual std::vector<double> computeGhostSolution(
    const std::vector<double> & U, const double & A) const override;
};

#endif /* BCEuler2PhaseGhostFree_H */
