#ifndef BCEuler2PhaseGhostDensityVelocity_H
#define BCEuler2PhaseGhostDensityVelocity_H

#include "BCEuler2PhaseGhost.h"

class BCEuler2PhaseGhostDensityVelocity : public BCEuler2PhaseGhost
{
public:
  BCEuler2PhaseGhostDensityVelocity(
    bool is_left, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap, const FluxEuler2Phase & flux,
    const double & a_liq,
    const double & r_liq, const double & u_liq,
    const double & r_vap, const double & u_vap);

protected:
  virtual std::vector<double> computeGhostSolution(
    const std::vector<double> & U, const double & A) const override;

  const double _a_liq;
  const double _r_liq;
  const double _u_liq;
  const double _r_vap;
  const double _u_vap;
};

#endif /* BCEuler2PhaseGhostDensityVelocity_H */
