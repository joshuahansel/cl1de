#ifndef ReconstructorEuler2PhaseSlopePUT_H
#define ReconstructorEuler2PhaseSlopePUT_H

#include "ReconstructorEuler2PhaseSlope.h"

class ReconstructorEuler2PhaseSlopePUT : public ReconstructorEuler2PhaseSlope
{
public:
  ReconstructorEuler2PhaseSlopePUT(
    const DoFHandlerEuler2Phase & dof_handler, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap);

protected:
  virtual void computeSlopeVariables(
    const double & aA_liq,
    const double & arA_liq,
    const double & aruA_liq,
    const double & arEA_liq,
    const double & arA_vap,
    const double & aruA_vap,
    const double & arEA_vap,
    const double & A,
    double & w1,
    double & w2,
    double & w3,
    double & w4,
    double & w5,
    double & w6,
    double & w7) const override;

  virtual void computeConservativeVariables(
    const double & w1,
    const double & w2,
    const double & w3,
    const double & w4,
    const double & w5,
    const double & w6,
    const double & w7,
    const double & A,
    double & aA_liq,
    double & arA_liq,
    double & aruA_liq,
    double & arEA_liq,
    double & arA_vap,
    double & aruA_vap,
    double & arEA_vap) const override;
};

#endif /* ReconstructorEuler2PhaseSlopePUT_H */
