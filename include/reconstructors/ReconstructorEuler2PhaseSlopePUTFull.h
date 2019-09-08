#ifndef ReconstructorEuler2PhaseSlopePUTFull_H
#define ReconstructorEuler2PhaseSlopePUTFull_H

#include "ReconstructorEuler2PhaseSlopePUT.h"

class EOS1Phase;

class ReconstructorEuler2PhaseSlopePUTFull : public ReconstructorEuler2PhaseSlopePUT
{
public:
  ReconstructorEuler2PhaseSlopePUTFull(
    const DoFHandlerEuler2Phase & dof_handler, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap);

protected:
  virtual double computeSlope(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R) const override;
};

#endif /* ReconstructorEuler2PhaseSlopePUTFull_H */
