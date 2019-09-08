#include "FunctionConstant.h"
#include "FunctionRiemann.h"
#include "EOS1PhaseIdealGas.h"
#include "BCEuler2PhaseGhostWall.h"
#include "ExecutionerEuler2Phase.h"
#include "FluxEuler2Phase.h"
#include "ICsEuler2PhaseRUP.h"
#include "ProblemEuler2Phase.h"
#include "RunParametersEuler2Phase.h"

int main(int argc, char* argv[])
{
  const double x_interface = 0.5;

  EOS1PhaseIdealGas eos_liq(1.4, 28.97e-3);
  EOS1PhaseIdealGas eos_vap(1.4, 28.97e-3);

  RunParametersEuler2Phase run_params(argc, argv, eos_liq, eos_vap);

  const FluxEuler2Phase & flux = run_params.getFlux();
  BCEuler2PhaseGhostWall left_bc(true, eos_liq, eos_vap, flux);
  BCEuler2PhaseGhostWall right_bc(false, eos_liq, eos_vap, flux);

  FunctionConstant ic_a(0.5);
  FunctionRiemann ic_r_liq(1.0, 0.125, x_interface);
  FunctionRiemann ic_r_vap(0.125, 1.0, x_interface);
  FunctionConstant ic_u(0.0);
  FunctionRiemann ic_p_liq(1.0, 0.1, x_interface);
  FunctionRiemann ic_p_vap(0.1, 1.0, x_interface);
  ICsEuler2PhaseRUP ics(ic_a, ic_r_liq, ic_u, ic_p_liq, ic_r_vap, ic_u, ic_p_vap, eos_liq, eos_vap);

  FunctionConstant area(1.0);
  ProblemEuler2Phase problem;
  problem.setAreaFunction(area);
  problem.setICs(ics);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);
  problem.setEOS(eos_liq, eos_vap);
  problem.setDomain(0.0, 1.0);
  problem.setDefaultEndTime(0.2);

  ExecutionerEuler2Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
