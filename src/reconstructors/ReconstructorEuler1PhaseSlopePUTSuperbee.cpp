#include "ReconstructorEuler1PhaseSlopePUTSuperbee.h"
#include "UtilsNumerics.h"

ReconstructorEuler1PhaseSlopePUTSuperbee::ReconstructorEuler1PhaseSlopePUTSuperbee(
  const DoFHandlerEuler1Phase & dof_handler, const EOS1Phase & eos)
  : ReconstructorEuler1PhaseSlopePUT(dof_handler, eos)
{
}

double ReconstructorEuler1PhaseSlopePUTSuperbee::computeSlope(
  const double & U_L, const double & U, const double & U_R,
  const double & x_L, const double & x, const double & x_R) const
{
  return Numerics::computeSlopeSuperbee(U_L, U, U_R, x_L, x, x_R);
}
