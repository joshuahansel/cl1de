#include "ReconstructorEuler1PhaseSlopeConservative.h"
#include "EOS1Phase.h"

ReconstructorEuler1PhaseSlopeConservative::ReconstructorEuler1PhaseSlopeConservative(
  const DoFHandlerEuler1Phase & dof_handler, const EOS1Phase & eos)
  : ReconstructorEuler1PhaseSlope(dof_handler, eos)
{
}

void ReconstructorEuler1PhaseSlopeConservative::computeSlopeVariables(
  const double & rA,
  const double & ruA,
  const double & rEA,
  const double & /*A*/,
  double & w_rA,
  double & w_ruA,
  double & w_rEA) const
{
  w_rA = rA;
  w_ruA = ruA;
  w_rEA = rEA;
}

void ReconstructorEuler1PhaseSlopeConservative::computeConservativeVariables(
  const double & w_rA,
  const double & w_ruA,
  const double & w_rEA,
  const double & /*A*/,
  double & rA,
  double & ruA,
  double & rEA) const
{
  rA = w_rA;
  ruA = w_ruA;
  rEA = w_rEA;
}
