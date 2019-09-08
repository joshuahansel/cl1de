#include "ReconstructorEuler2PhaseSlopePUTMC.h"
#include "UtilsNumerics.h"

ReconstructorEuler2PhaseSlopePUTMC::ReconstructorEuler2PhaseSlopePUTMC(
  const DoFHandlerEuler2Phase & dof_handler, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap)
  : ReconstructorEuler2PhaseSlopePUT(dof_handler, eos_liq, eos_vap)
{
}

double ReconstructorEuler2PhaseSlopePUTMC::computeSlope(
  const double & U_L, const double & U, const double & U_R,
  const double & x_L, const double & x, const double & x_R) const
{
  return Numerics::computeSlopeMC(U_L, U, U_R, x_L, x, x_R);
}
