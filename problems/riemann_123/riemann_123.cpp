#include "BCEuler1PhaseFree.h"
#include "EOS1PhaseIdealGas.h"
#include "ExecutionerEuler1Phase.h"
#include "FunctionConstant.h"
#include "FunctionRiemann.h"
#include "ICsEuler1PhaseRUP.h"
#include "ProblemEuler1Phase.h"
#include "RunParametersEuler1Phase.h"

int main(int argc, char* argv[])
{
  ProblemEuler1Phase problem;
  problem.setDomain(0.0, 1.0);
  problem.setDefaultEndTime(0.15);

  EOS1PhaseIdealGas eos(1.4, 28.97e-3);
  problem.setEOS(eos);

  RunParametersEuler1Phase run_params(argc, argv, eos);

  const double x_interface = 0.5;
  FunctionConstant area(1.0);
  FunctionRiemann velocity_ic(-2.0, 2.0, x_interface);
  FunctionConstant pressure_ic(0.4);
  FunctionConstant density_ic(1.0);
  ICsEuler1PhaseRUP ics(density_ic, velocity_ic, pressure_ic, eos);
  problem.setAreaFunction(area);
  problem.setICs(ics);

  BCEuler1PhaseFree left_bc(true, eos);
  BCEuler1PhaseFree right_bc(false, eos);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);

  ExecutionerEuler1Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
