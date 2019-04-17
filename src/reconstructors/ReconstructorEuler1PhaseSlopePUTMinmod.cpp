#include "ReconstructorEuler1PhaseSlopePUTMinmod.h"
#include "UtilsNumerics.h"

ReconstructorEuler1PhaseSlopePUTMinmod::ReconstructorEuler1PhaseSlopePUTMinmod(const EOS1Phase & eos)
  : ReconstructorEuler1PhaseSlopePUT(eos)
{
}

double ReconstructorEuler1PhaseSlopePUTMinmod::computeSlope(
  const double & U_L, const double & U, const double & U_R,
  const double & x_L, const double & x, const double & x_R) const
{
  return Numerics::computeSlopeMinmod(U_L, U, U_R, x_L, x, x_R);
}
