#ifndef ReconstructorEuler1PhaseSlopeRUP_H
#define ReconstructorEuler1PhaseSlopeRUP_H

#include "ReconstructorEuler1PhaseSlope.h"

class ReconstructorEuler1PhaseSlopeRUP : public ReconstructorEuler1PhaseSlope
{
public:
  ReconstructorEuler1PhaseSlopeRUP(const EOS1Phase & eos);

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

#endif /* ReconstructorEuler1PhaseSlopeRUP_H */
