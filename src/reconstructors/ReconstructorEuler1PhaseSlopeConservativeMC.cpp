#include "ReconstructorEuler1PhaseSlopeConservativeMC.h"
#include "UtilsNumerics.h"

ReconstructorEuler1PhaseSlopeConservativeMC::ReconstructorEuler1PhaseSlopeConservativeMC(const EOS1Phase & eos)
  : ReconstructorEuler1PhaseSlopeConservative(eos)
{
}

double ReconstructorEuler1PhaseSlopeConservativeMC::computeSlope(
  const double & U_L, const double & U, const double & U_R,
  const double & x_L, const double & x, const double & x_R) const
{
  return Numerics::computeSlopeMC(U_L, U, U_R, x_L, x, x_R);
}
