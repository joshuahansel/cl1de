#ifndef ReconstructorEuler2PhaseSlopePUTSuperbee_H
#define ReconstructorEuler2PhaseSlopePUTSuperbee_H

#include "ReconstructorEuler2PhaseSlopePUT.h"

class EOS1Phase;

class ReconstructorEuler2PhaseSlopePUTSuperbee : public ReconstructorEuler2PhaseSlopePUT
{
public:
  ReconstructorEuler2PhaseSlopePUTSuperbee(
    const DoFHandlerEuler2Phase & dof_handler, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap);

protected:
  virtual double computeSlope(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R) const override;
};

#endif /* ReconstructorEuler2PhaseSlopePUTSuperbee_H */
