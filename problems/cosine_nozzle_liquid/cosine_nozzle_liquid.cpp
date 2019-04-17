#include "BCEuler1PhaseGhostPressure.h"
#include "BCEuler1PhaseGhostStagnation.h"
#include "EOS1PhaseStiffenedGas.h"
#include "ExecutionerEuler1Phase.h"
#include "FluxEuler1Phase.h"
#include "FunctionConstant.h"
#include "FunctionCosineHump.h"
#include "FunctionDensityFromPressureTemperature.h"
#include "ProblemEuler1Phase.h"
#include "RunParametersEuler1Phase.h"

int main(int argc, char* argv[])
{
  ProblemEuler1Phase problem;
  problem.setDomain(0.0, 1.0);
  problem.setDefaultEndTime(1.0);

  EOS1PhaseStiffenedGas eos(2.35, 1816, 1e9, -1.167e6, 0);
  problem.setEOS(eos);

  RunParametersEuler1Phase run_params(argc, argv, eos);

  const double p0_inlet = 1.0e6;
  const double T0_inlet = 1073;
  const double p_outlet = 0.5e6;
  FunctionCosineHump area(1.0, 0.5, 0.5, 1.0);
  FunctionConstant velocity_ic(0.0);
  FunctionConstant pressure_ic(p0_inlet);
  FunctionConstant temperature_ic(T0_inlet);
  FunctionDensityFromPressureTemperature density_ic(eos, pressure_ic, temperature_ic);
  problem.setAreaFunction(area);
  problem.setICDensity(density_ic);
  problem.setICVelocity(velocity_ic);
  problem.setICPressure(pressure_ic);

  const FluxEuler1Phase & flux = run_params.getFlux();
  BCEuler1PhaseGhostStagnation left_bc(true, eos, flux, p0_inlet, T0_inlet);
  BCEuler1PhaseGhostPressure right_bc(false, eos, flux, p_outlet);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);

  ExecutionerEuler1Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
