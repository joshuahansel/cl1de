#ifndef BCEuler2PhaseGhostSetRUP_H
#define BCEuler2PhaseGhostSetRUP_H

#include "BCEuler2PhaseGhost.h"

class BCEuler2PhaseGhostSetRUP : public BCEuler2PhaseGhost
{
public:
  BCEuler2PhaseGhostSetRUP(
    bool is_left, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap, const FluxEuler2Phase & flux,
    const double & a_liq,
    const double & r_liq, const double & u_liq, const double & p_liq,
    const double & r_vap, const double & u_vap, const double & p_vap);

protected:
  virtual std::vector<double> computeGhostSolution(
    const std::vector<double> & U, const double & A) const override;

  const double _a_liq;
  const double _r_liq;
  const double _u_liq;
  const double _p_liq;
  const double _r_vap;
  const double _u_vap;
  const double _p_vap;
};

#endif /* BCEuler2PhaseGhostSetRUP_H */
