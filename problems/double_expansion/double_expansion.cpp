#include "FunctionConstant.h"
#include "FunctionRiemann.h"
#include "BCEuler2PhaseGhostFree.h"
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

  RunParametersEuler2Phase run_params(argc, argv, eos_liq, eos_vap);

  const double a_liq = 0.5;
  const double r_liq = 1000.0;
  const double r_vap = 10.0;
  const double u_L = -1000.0;
  const double u_R = 1000.0;
  const double p = 0.1e6;

  const FluxEuler2Phase & flux = run_params.getFlux();
  BCEuler2PhaseGhostFree left_bc(true, eos_liq, eos_vap, flux);
  BCEuler2PhaseGhostFree right_bc(false, eos_liq, eos_vap, flux);

  const double x_interface = 0.5;
  FunctionConstant ic_a_liq(a_liq);
  FunctionConstant ic_r_liq(r_liq);
  FunctionConstant ic_r_vap(r_vap);
  FunctionRiemann ic_u(u_L, u_R, x_interface);
  FunctionConstant ic_p(p);
  ICsEuler2PhaseRUP ics(ic_a_liq, ic_r_liq, ic_u, ic_p, ic_r_vap, ic_u, ic_p, eos_liq, eos_vap);

  FunctionConstant area(1.0);
  ProblemEuler2Phase problem;
  problem.setAreaFunction(area);
  problem.setICs(ics);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);
  problem.setEOS(eos_liq, eos_vap);
  problem.setDomain(0.0, 1.0);
  problem.setDefaultEndTime(207e-6);

  ExecutionerEuler2Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
