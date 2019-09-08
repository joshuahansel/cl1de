#ifndef BCEuler1PhaseGhostWall_H
#define BCEuler1PhaseGhostWall_H

#include "BCEuler1PhaseGhost.h"

class BCEuler1PhaseGhostWall : public BCEuler1PhaseGhost
{
public:
  BCEuler1PhaseGhostWall(
    bool is_left, const EOS1Phase & eos, const FluxEuler1Phase & flux);

protected:
  virtual std::vector<double> computeGhostSolution(
    const std::vector<double> & U, const double & A) const override;
};

#endif /* BCEuler1PhaseGhostWall_H */
