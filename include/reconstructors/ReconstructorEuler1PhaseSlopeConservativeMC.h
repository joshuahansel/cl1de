#ifndef ReconstructorEuler1PhaseSlopeConservativeMC_H
#define ReconstructorEuler1PhaseSlopeConservativeMC_H

#include "ReconstructorEuler1PhaseSlopeConservative.h"

class EOS1Phase;

class ReconstructorEuler1PhaseSlopeConservativeMC : public ReconstructorEuler1PhaseSlopeConservative
{
public:
  ReconstructorEuler1PhaseSlopeConservativeMC(const EOS1Phase & eos);

protected:
  virtual double computeSlope(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R) const override;
};

#endif /* ReconstructorEuler1PhaseSlopeConservativeMC_H */
