#ifndef BCEuler2PhaseGhost_H
#define BCEuler2PhaseGhost_H

#include "BCEuler2Phase.h"

class EOS1Phase;
class FluxEuler2Phase;

class BCEuler2PhaseGhost : public BCEuler2Phase
{
public:
  BCEuler2PhaseGhost(
    bool is_left, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap, const FluxEuler2Phase & flux);

  virtual std::vector<double> computeFlux(
    const std::vector<double> & U, const double & A) const override;

protected:
  virtual std::vector<double> computeGhostSolution(
    const std::vector<double> & U, const double & A) const = 0;

  const EOS1Phase & _eos_liq;
  const EOS1Phase & _eos_vap;
  const FluxEuler2Phase & _flux;
};

#endif /* BCEuler2PhaseGhost_H */
