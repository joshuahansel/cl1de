#include "ReconstructorEuler1PhaseSlopeRUP.h"
#include "EOS1Phase.h"

ReconstructorEuler1PhaseSlopeRUP::ReconstructorEuler1PhaseSlopeRUP(
  const DoFHandlerEuler1Phase & dof_handler, const EOS1Phase & eos)
  : ReconstructorEuler1PhaseSlope(dof_handler, eos)
{
}

void ReconstructorEuler1PhaseSlopeRUP::computeSlopeVariables(
  const double & rA,
  const double & ruA,
  const double & rEA,
  const double & A,
  double & r,
  double & u,
  double & p) const
{
  r = rA / A;
  u = ruA / rA;
  const double E = rEA / rA;
  const double e = E - 0.5 * u * u;
  p = _eos.p_from_r_e(r, e);
}

void ReconstructorEuler1PhaseSlopeRUP::computeConservativeVariables(
  const double & r,
  const double & u,
  const double & p,
  const double & A,
  double & rA,
  double & ruA,
  double & rEA) const
{
  const double e = _eos.e_from_p_r(p, r);
  const double E = e + 0.5 * u * u;

  rA = r * A;
  ruA = rA * u;
  rEA = rA * E;
}
