#include "FunctionConstant.h"
#include "FunctionRiemann.h"
#include "BCEuler1PhaseFree.h"
#include "EOS1PhaseIdealGas.h"
#include "ProblemEuler1Phase.h"
#include "RunParametersEuler1Phase.h"
#include "ExecutionerEuler1Phase.h"

int main(int argc, char* argv[])
{
  const double x_interface = 0.5;
  FunctionConstant area(1.0);
  FunctionRiemann density_ic(1.0, 0.125, x_interface);
  FunctionConstant velocity_ic(0.0);
  FunctionRiemann pressure_ic(1.0, 0.1, x_interface);
  EOS1PhaseIdealGas eos(1.4, 28.97e-3);
  BCEuler1PhaseFree left_bc(true, eos);
  BCEuler1PhaseFree right_bc(false, eos);

  ProblemEuler1Phase problem;
  problem.setAreaFunction(area);
  problem.setICDensity(density_ic);
  problem.setICVelocity(velocity_ic);
  problem.setICPressure(pressure_ic);
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
