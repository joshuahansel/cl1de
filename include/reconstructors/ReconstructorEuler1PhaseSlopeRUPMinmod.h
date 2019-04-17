#ifndef ReconstructorEuler1PhaseSlopeRUPMinmod_H
#define ReconstructorEuler1PhaseSlopeRUPMinmod_H

#include "ReconstructorEuler1PhaseSlopeRUP.h"

class EOS1Phase;

class ReconstructorEuler1PhaseSlopeRUPMinmod : public ReconstructorEuler1PhaseSlopeRUP
{
public:
  ReconstructorEuler1PhaseSlopeRUPMinmod(const EOS1Phase & eos);

protected:
  virtual double computeSlope(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R) const override;
};

#endif /* ReconstructorEuler1PhaseSlopeRUPMinmod_H */
