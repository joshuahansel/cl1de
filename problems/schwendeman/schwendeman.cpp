#include "FunctionConstant.h"
#include "FunctionRiemann.h"
#include "BCEuler2PhaseGhostFree.h"
#include "EOS1PhaseIdealGas.h"
#include "EOS1PhaseStiffenedGas.h"
#include "ExecutionerEuler2Phase.h"
#include "FluxEuler2Phase.h"
#include "ICsEuler2PhaseRUP.h"
#include "ProblemEuler2Phase.h"
#include "RunParametersEuler2Phase.h"

int main(int argc, char* argv[])
{
  EOS1PhaseStiffenedGas eos_liq(1.4, 1.0, 0, 0, 0);
  EOS1PhaseIdealGas eos_vap(1.4, 1.0);

  RunParametersEuler2Phase run_params(argc, argv, eos_liq, eos_vap);

  const double a_liq_L = 0.8;
  const double a_liq_R = 0.3;
  const double r_liq = 1.0;
  const double r_vap_L = 0.2;
  const double r_vap_R = 1.0;
  const double u = 0.0;
  const double p_liq = 1.0;
  const double p_vap_L = 0.3;
  const double p_vap_R = 1.0;

  const FluxEuler2Phase & flux = run_params.getFlux();
  BCEuler2PhaseGhostFree left_bc(true, eos_liq, eos_vap, flux);
  BCEuler2PhaseGhostFree right_bc(false, eos_liq, eos_vap, flux);

  const double x_interface = 0.0;
  FunctionRiemann ic_a_liq(a_liq_L, a_liq_R, x_interface);
  FunctionConstant ic_r_liq(r_liq);
  FunctionRiemann ic_r_vap(r_vap_L, r_vap_R, x_interface);
  FunctionConstant ic_u(u);
  FunctionConstant ic_p_liq(p_liq);
  FunctionRiemann ic_p_vap(p_vap_L, p_vap_R, x_interface);
  ICsEuler2PhaseRUP ics(ic_a_liq, ic_r_liq, ic_u, ic_p_liq, ic_r_vap, ic_u, ic_p_vap, eos_liq, eos_vap);

  FunctionConstant area(1.0);
  ProblemEuler2Phase problem;
  problem.setAreaFunction(area);
  problem.setICs(ics);
  problem.setBCLeft(left_bc);
  problem.setBCRight(right_bc);
  problem.setEOS(eos_liq, eos_vap);
  problem.setDomain(-0.5, 0.5);
  problem.setDefaultEndTime(0.2);

  ExecutionerEuler2Phase executioner(problem, run_params);
  executioner.run();

  return 0;
}
