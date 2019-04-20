#include "ReconstructorEuler1PhaseSlopePUT.h"
#include "EOS1Phase.h"

ReconstructorEuler1PhaseSlopePUT::ReconstructorEuler1PhaseSlopePUT(
  const DoFHandlerEuler1Phase & dof_handler, const EOS1Phase & eos)
  : ReconstructorEuler1PhaseSlope(dof_handler, eos)
{
}

void ReconstructorEuler1PhaseSlopePUT::computeSlopeVariables(
  const double & rA,
  const double & ruA,
  const double & rEA,
  const double & A,
  double & p,
  double & u,
  double & T) const
{
  const double r = rA / A;
  u = ruA / rA;
  const double E = rEA / rA;
  const double e = E - 0.5 * u * u;
  p = _eos.p_from_r_e(r, e);
  T = _eos.T_from_r_e(r, e);
}

void ReconstructorEuler1PhaseSlopePUT::computeConservativeVariables(
  const double & p,
  const double & u,
  const double & T,
  const double & A,
  double & rA,
  double & ruA,
  double & rEA) const
{
  const double r = _eos.r_from_p_T(p, T);
  const double e = _eos.e_from_p_r(p, r);
  const double E = e + 0.5 * u * u;

  rA = r * A;
  ruA = rA * u;
  rEA = rA * E;
}
