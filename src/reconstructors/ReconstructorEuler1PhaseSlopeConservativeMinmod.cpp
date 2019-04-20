#include "ReconstructorEuler1PhaseSlopeConservativeMinmod.h"
#include "UtilsNumerics.h"

ReconstructorEuler1PhaseSlopeConservativeMinmod::ReconstructorEuler1PhaseSlopeConservativeMinmod(
  const DoFHandlerEuler1Phase & dof_handler, const EOS1Phase & eos)
  : ReconstructorEuler1PhaseSlopeConservative(dof_handler, eos)
{
}

double ReconstructorEuler1PhaseSlopeConservativeMinmod::computeSlope(
  const double & U_L, const double & U, const double & U_R,
  const double & x_L, const double & x, const double & x_R) const
{
  return Numerics::computeSlopeMinmod(U_L, U, U_R, x_L, x, x_R);
}
