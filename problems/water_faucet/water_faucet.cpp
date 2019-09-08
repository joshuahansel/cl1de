#include "FunctionConstant.h"
#include "BCEuler2PhaseGhostDensityVelocity.h"
#include "BCEuler2PhaseGhostPressure.h"
#include "EOS1PhaseStiffenedGas.h"
#include "ExecutionerEuler2Phase.h"
#include "FluxEuler2Phase.h"
#include "ICsEuler2PhaseRUP.h"
#include "ProblemEuler2Phase.h"
#include "RunParametersEuler2Phase.h"

int main(int argc, char* argv[])
{
  EOS1PhaseStiffenedGas eos_liq(4.4, 1.0, 6.0e8, 0, 0);
  EOS1PhaseStiffenedGas eos_vap(1.4, 1.0, 0, 0, 0);

  const double a_liq = 0.8;
  const double r_liq = 1000.0;
  const double r_vap = 1.0;
  const double u_liq = 10.0;
  const double u_vap = 0.0;
  const double p = 1e5;

  RunParametersEuler2Phase run_params(argc, argv, eos_liq, eos_vap);

  const FluxEuler2Phase & flux = run_params.getFlux();
  BCEuler2PhaseGhostDensityVelocity left_bc(true, eos_liq, eos_vap, flux, a_liq, r_liq, u_liq, r_vap, u_vap);
  BCEuler2PhaseGhostPressure right_bc(false, eos_liq, eos_vap, flux, p);

  FunctionConstant ic_a_liq(a_liq);
  FunctionConstant ic_r_liq(r_liq);
  FunctionConstant ic_r_vap(r_vap);
  FunctionConstant ic_u_liq(u_liq);
  FunctionConstant ic_u_vap(u_vap);
  FunctionConstant ic_p(p);
  ICsEuler2PhaseRUP ics(ic_a_liq, ic_r_liq, ic_u_liq, ic_p, ic_r_vap, ic_u_vap, ic_p, eos_liq, eos_vap);

  FunctionConstant area(1.0);
  ProblemEuler2Phase problem;
  problem.setAreaFunction(area);
  problem.setICs(ics);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);
  problem.setEOS(eos_liq, eos_vap);
  problem.setGravity(10.0);
  problem.setDomain(0.0, 12.0);
  problem.setDefaultEndTime(0.4);

  ExecutionerEuler2Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
