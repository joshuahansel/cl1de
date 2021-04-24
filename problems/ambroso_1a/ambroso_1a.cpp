#include "FunctionConstant.h"
#include "FunctionRiemann.h"
#include "BCEuler2PhaseGhostFree.h"
#include "EOS1PhaseIdealGas.h"
#include "ExecutionerEuler2Phase.h"
#include "FluxEuler2Phase.h"
#include "ICsEuler2PhaseRUP.h"
#include "ProblemEuler2Phase.h"
#include "RunParametersEuler2Phase.h"

int main(int argc, char* argv[])
{
  EOS1PhaseIdealGas eos_liq(1.4, 1.0);
  EOS1PhaseIdealGas eos_vap(1.4, 1.0);

  RunParametersEuler2Phase run_params(argc, argv, eos_liq, eos_vap);

  const double a_liq_L = 0.9;
  const double a_liq_R = 0.5;
  const double r_liq_L = 1.0;
  const double r_liq_R = 0.125;
  const double r_vap_L = 1.0;
  const double r_vap_R = 0.125;
  const double u_liq_L = 100.0;
  const double u_liq_R = 100.0;
  const double u_vap_L = 100.0;
  const double u_vap_R = 100.0;
  const double p_liq_L = 1e5;
  const double p_liq_R = 1e5;
  const double p_vap_L = 1e5;
  const double p_vap_R = 1e5;

  const FluxEuler2Phase & flux = run_params.getFlux();
  BCEuler2PhaseGhostFree left_bc(true, eos_liq, eos_vap, flux);
  BCEuler2PhaseGhostFree right_bc(false, eos_liq, eos_vap, flux);

  const double x_interface = 0.5;
  FunctionRiemann ic_a_liq(a_liq_L, a_liq_R, x_interface);
  FunctionRiemann ic_r_liq(r_liq_L, r_liq_R, x_interface);
  FunctionRiemann ic_r_vap(r_vap_L, r_vap_R, x_interface);
  FunctionRiemann ic_u_liq(u_liq_L, u_liq_R, x_interface);
  FunctionRiemann ic_u_vap(u_vap_L, u_vap_R, x_interface);
  FunctionRiemann ic_p_liq(p_liq_L, p_liq_R, x_interface);
  FunctionRiemann ic_p_vap(p_vap_L, p_vap_R, x_interface);
  ICsEuler2PhaseRUP ics(ic_a_liq, ic_r_liq, ic_u_liq, ic_p_liq, ic_r_vap, ic_u_vap, ic_p_vap, eos_liq, eos_vap);

  FunctionConstant area(1.0);
  ProblemEuler2Phase problem;
  problem.setAreaFunction(area);
  problem.setICs(ics);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);
  problem.setEOS(eos_liq, eos_vap);
  problem.setDomain(0.0, 1.0);
  problem.setDefaultEndTime(3.0);

  ExecutionerEuler2Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
