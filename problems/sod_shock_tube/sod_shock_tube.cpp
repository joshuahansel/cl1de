#include "FunctionConstant.h"
#include "FunctionRiemann.h"
#include "BCEuler1PhaseFree.h"
#include "EOS1PhaseIdealGas.h"
#include "ICsEuler1PhaseRUP.h"
#include "ProblemEuler1Phase.h"
#include "RunParametersEuler1Phase.h"
#include "ExecutionerEuler1Phase.h"

int main(int argc, char* argv[])
{
  const double x_interface = 0.5;

  EOS1PhaseIdealGas eos(1.4, 28.97e-3);
  FunctionConstant area(1.0);
  FunctionRiemann density_ic(1.0, 0.125, x_interface);
  FunctionConstant velocity_ic(0.0);
  FunctionRiemann pressure_ic(1.0, 0.1, x_interface);
  ICsEuler1PhaseRUP ics(density_ic, velocity_ic, pressure_ic, eos);
  BCEuler1PhaseFree left_bc(true, eos);
  BCEuler1PhaseFree right_bc(false, eos);

  ProblemEuler1Phase problem;
  problem.setAreaFunction(area);
  problem.setICs(ics);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);
  problem.setEOS(eos);
  problem.setDomain(0.0, 1.0);
  problem.setDefaultEndTime(0.2);

  RunParametersEuler1Phase run_params(argc, argv, eos);
  ExecutionerEuler1Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
