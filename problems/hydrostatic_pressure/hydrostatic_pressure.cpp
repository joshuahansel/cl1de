#include "FunctionConstant.h"
#include "BCEuler1PhaseGhostWall.h"
#include "EOS1PhaseStiffenedGas.h"
#include "ExecutionerEuler1Phase.h"
#include "FluxEuler1Phase.h"
#include "ICsEuler1PhaseRUP.h"
#include "ProblemEuler1Phase.h"
#include "RunParametersEuler1Phase.h"

int main(int argc, char* argv[])
{
  EOS1PhaseStiffenedGas eos(1.4, 1.0, 0, 0, 0);

  RunParametersEuler1Phase run_params(argc, argv, eos);

  const FluxEuler1Phase & flux = run_params.getFlux();
  BCEuler1PhaseGhostWall left_bc(true, eos, flux);
  BCEuler1PhaseGhostWall right_bc(false, eos, flux);

  FunctionConstant ic_r(1.0);
  FunctionConstant ic_u(0.0);
  FunctionConstant ic_p(1.0e5);
  ICsEuler1PhaseRUP ics(ic_r, ic_u, ic_p, eos);

  FunctionConstant area(1.0);
  ProblemEuler1Phase problem;
  problem.setAreaFunction(area);
  problem.setICs(ics);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);
  problem.setEOS(eos);
  problem.setGravity(-10.0);
  problem.setDomain(0.0, 1.0);
  problem.setDefaultEndTime(1.6);

  ExecutionerEuler1Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
