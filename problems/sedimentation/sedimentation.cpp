#include "FunctionConstant.h"
#include "BCEuler2PhaseGhostWall.h"
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

  const FluxEuler2Phase & flux = run_params.getFlux();
  BCEuler2PhaseGhostWall left_bc(true, eos_liq, eos_vap, flux);
  BCEuler2PhaseGhostWall right_bc(false, eos_liq, eos_vap, flux);

  FunctionConstant ic_a_liq(0.5);
  FunctionConstant ic_r_liq(1000.0);
  FunctionConstant ic_r_vap(1.0);
  FunctionConstant ic_u(0.0);
  FunctionConstant ic_p(1.0e5);
  ICsEuler2PhaseRUP ics(ic_a_liq, ic_r_liq, ic_u, ic_p, ic_r_vap, ic_u, ic_p, eos_liq, eos_vap);

  FunctionConstant area(1.0);
  ProblemEuler2Phase problem;
  problem.setAreaFunction(area);
  problem.setICs(ics);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);
  problem.setEOS(eos_liq, eos_vap);
  problem.setGravity(-10.0);
  problem.setDomain(0.0, 1.0);
  problem.setDefaultEndTime(1.6);

  ExecutionerEuler2Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
