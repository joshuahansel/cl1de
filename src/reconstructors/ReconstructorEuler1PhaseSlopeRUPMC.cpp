#include "ReconstructorEuler1PhaseSlopeRUPMC.h"
#include "UtilsNumerics.h"

ReconstructorEuler1PhaseSlopeRUPMC::ReconstructorEuler1PhaseSlopeRUPMC(const EOS1Phase & eos)
  : ReconstructorEuler1PhaseSlopeRUP(eos)
{
}

double ReconstructorEuler1PhaseSlopeRUPMC::computeSlope(
  const double & U_L, const double & U, const double & U_R,
  const double & x_L, const double & x, const double & x_R) const
{
  return Numerics::computeSlopeMC(U_L, U, U_R, x_L, x, x_R);
}
