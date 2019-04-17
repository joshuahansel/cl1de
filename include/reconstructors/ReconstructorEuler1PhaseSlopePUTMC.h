#ifndef ReconstructorEuler1PhaseSlopePUTMC_H
#define ReconstructorEuler1PhaseSlopePUTMC_H

#include "ReconstructorEuler1PhaseSlopePUT.h"

class EOS1Phase;

class ReconstructorEuler1PhaseSlopePUTMC : public ReconstructorEuler1PhaseSlopePUT
{
public:
  ReconstructorEuler1PhaseSlopePUTMC(const EOS1Phase & eos);

protected:
  virtual double computeSlope(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R) const override;
};

#endif /* ReconstructorEuler1PhaseSlopePUTMC_H */
