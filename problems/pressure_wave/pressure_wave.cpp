#include "BCEuler1PhaseGhostPressure.h"
#include "EOS1PhaseIdealGas.h"
#include "ExecutionerEuler1Phase.h"
#include "FluxEuler1Phase.h"
#include "FunctionConstant.h"
#include "FunctionDensityFromPressureTemperature.h"
#include "ProblemEuler1Phase.h"
#include "RunParametersEuler1Phase.h"

int main(int argc, char* argv[])
{
  ProblemEuler1Phase problem;
  problem.setDomain(0.0, 4.0);
  problem.setDefaultEndTime(0.006);

  EOS1PhaseIdealGas eos(1.4, 28.97e-3);
  problem.setEOS(eos);

  RunParametersEuler1Phase run_params(argc, argv, eos);

  const double p_inlet = 200.0e3;
  const double p_initial = 100.0e3;
  FunctionConstant area(1.0);
  FunctionConstant velocity_ic(0.0);
  FunctionConstant pressure_ic(p_initial);
  FunctionConstant temperature_ic(300.0);
  FunctionDensityFromPressureTemperature density_ic(eos, pressure_ic, temperature_ic);
  problem.setAreaFunction(area);
  problem.setICDensity(density_ic);
  problem.setICVelocity(velocity_ic);
  problem.setICPressure(pressure_ic);

  const FluxEuler1Phase & flux = run_params.getFlux();
  BCEuler1PhaseGhostPressure left_bc(true, eos, flux, p_inlet);
  BCEuler1PhaseGhostPressure right_bc(false, eos, flux, p_initial);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);

  ExecutionerEuler1Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
