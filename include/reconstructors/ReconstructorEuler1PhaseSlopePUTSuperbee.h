#ifndef ReconstructorEuler1PhaseSlopePUTSuperbee_H
#define ReconstructorEuler1PhaseSlopePUTSuperbee_H

#include "ReconstructorEuler1PhaseSlopePUT.h"

class EOS1Phase;

class ReconstructorEuler1PhaseSlopePUTSuperbee : public ReconstructorEuler1PhaseSlopePUT
{
public:
  ReconstructorEuler1PhaseSlopePUTSuperbee(const EOS1Phase & eos);

protected:
  virtual double computeSlope(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R) const override;
};

#endif /* ReconstructorEuler1PhaseSlopePUTSuperbee_H */
