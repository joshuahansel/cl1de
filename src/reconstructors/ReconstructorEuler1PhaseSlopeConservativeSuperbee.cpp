#include "ReconstructorEuler1PhaseSlopeConservativeSuperbee.h"
#include "UtilsNumerics.h"

ReconstructorEuler1PhaseSlopeConservativeSuperbee::ReconstructorEuler1PhaseSlopeConservativeSuperbee(
  const DoFHandlerEuler1Phase & dof_handler, const EOS1Phase & eos)
  : ReconstructorEuler1PhaseSlopeConservative(dof_handler, eos)
{
}

double ReconstructorEuler1PhaseSlopeConservativeSuperbee::computeSlope(
  const double & U_L, const double & U, const double & U_R,
  const double & x_L, const double & x, const double & x_R) const
{
  return Numerics::computeSlopeSuperbee(U_L, U, U_R, x_L, x, x_R);
}
