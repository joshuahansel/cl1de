#ifndef ReconstructorEuler1PhaseSlopeRUPSuperbee_H
#define ReconstructorEuler1PhaseSlopeRUPSuperbee_H

#include "ReconstructorEuler1PhaseSlopeRUP.h"

class EOS1Phase;

class ReconstructorEuler1PhaseSlopeRUPSuperbee : public ReconstructorEuler1PhaseSlopeRUP
{
public:
  ReconstructorEuler1PhaseSlopeRUPSuperbee(
    const DoFHandlerEuler1Phase & dof_handler, const EOS1Phase & eos);

protected:
  virtual double computeSlope(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R) const override;
};

#endif /* ReconstructorEuler1PhaseSlopeRUPSuperbee_H */
