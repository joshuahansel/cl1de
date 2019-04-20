#ifndef ReconstructorEuler1PhaseSlopePUTMinmod_H
#define ReconstructorEuler1PhaseSlopePUTMinmod_H

#include "ReconstructorEuler1PhaseSlopePUT.h"

class EOS1Phase;

class ReconstructorEuler1PhaseSlopePUTMinmod : public ReconstructorEuler1PhaseSlopePUT
{
public:
  ReconstructorEuler1PhaseSlopePUTMinmod(
    const DoFHandlerEuler1Phase & dof_handler, const EOS1Phase & eos);

protected:
  virtual double computeSlope(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R) const override;
};

#endif /* ReconstructorEuler1PhaseSlopePUTMinmod_H */
