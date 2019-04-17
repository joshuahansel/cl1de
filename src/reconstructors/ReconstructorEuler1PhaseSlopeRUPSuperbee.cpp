#include "ReconstructorEuler1PhaseSlopeRUPSuperbee.h"
#include "UtilsNumerics.h"

ReconstructorEuler1PhaseSlopeRUPSuperbee::ReconstructorEuler1PhaseSlopeRUPSuperbee(const EOS1Phase & eos)
  : ReconstructorEuler1PhaseSlopeRUP(eos)
{
}

double ReconstructorEuler1PhaseSlopeRUPSuperbee::computeSlope(
  const double & U_L, const double & U, const double & U_R,
  const double & x_L, const double & x, const double & x_R) const
{
  return Numerics::computeSlopeSuperbee(U_L, U, U_R, x_L, x, x_R);
}
