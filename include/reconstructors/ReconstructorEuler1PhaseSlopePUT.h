#ifndef ReconstructorEuler1PhaseSlopePUT_H
#define ReconstructorEuler1PhaseSlopePUT_H

#include "ReconstructorEuler1PhaseSlope.h"

class ReconstructorEuler1PhaseSlopePUT : public ReconstructorEuler1PhaseSlope
{
public:
  ReconstructorEuler1PhaseSlopePUT(
    const DoFHandlerEuler1Phase & dof_handler, const EOS1Phase & eos);

protected:
  virtual void computeSlopeVariables(
    const double & rA,
    const double & ruA,
    const double & rEA,
    const double & A,
    double & w1,
    double & w2,
    double & w3) const override;

  virtual void computeConservativeVariables(
    const double & w1,
    const double & w2,
    const double & w3,
    const double & A,
    double & rA,
    double & ruA,
    double & rEA) const override;
};

#endif /* ReconstructorEuler1PhaseSlopePUT_H */
