#include "ExecutionerEuler1Phase.h"
#include "ProblemEuler1Phase.h"
#include "RunParametersEuler1Phase.h"
#include "BCEuler1Phase.h"
#include "TimeStepSizer.h"
#include "ReconstructorEuler1Phase.h"
#include "FluxEuler1Phase.h"
#include "Function.h"
#include "EOS1Phase.h"
#include "CSVOutput.h"
#include "UtilsNumerics.h"

#include <iostream>
#include <vector>

ExecutionerEuler1Phase::ExecutionerEuler1Phase(
  const ProblemEuler1Phase & problem,
  const RunParametersEuler1Phase & run_params)
  : Executioner(problem, run_params),
    _problem(problem),
    _run_params(run_params),

    _eos(problem.getEOS()),
    _A_fn(problem.getAreaFunction()),
    _r_ic_fn(problem.getICDensity()),
    _u_ic_fn(problem.getICVelocity()),
    _p_ic_fn(problem.getICPressure()),

    _A_node(computeAreaNode()),
    _A_elem(computeAreaElem())
{
}

void ExecutionerEuler1Phase::run()
{
  std::vector<std::vector<double>> rA(_n_stages + 1, std::vector<double>(_n_elems, 0));
  std::vector<std::vector<double>> ruA(_n_stages + 1, std::vector<double>(_n_elems, 0));
  std::vector<std::vector<double>> rEA(_n_stages + 1, std::vector<double>(_n_elems, 0));

  // initialize
  double t = 0.0;

  for (unsigned int i = 0; i < _n_elems; i++)
  {
    const double r = _r_ic_fn.value(_x_elem[i], t);
    const double u = _u_ic_fn.value(_x_elem[i], t);
    const double p = _p_ic_fn.value(_x_elem[i], t);

    const double e = _eos.e_from_p_r(p, r);
    const double E = e + 0.5 * u * u;

    rA[0][i] = r * _A_elem[i];
    ruA[0][i] = r * u * _A_elem[i];
    rEA[0][i] = r * E * _A_elem[i];
  }

  const BCEuler1Phase & bc_left = _problem.getBCLeft();
  const BCEuler1Phase & bc_right = _problem.getBCRight();

  const ReconstructorEuler1Phase & reconstructor = _run_params.getReconstructor();
  const FluxEuler1Phase & flux = _run_params.getFlux();

  std::vector<double> rA_L(_n_elems);
  std::vector<double> ruA_L(_n_elems);
  std::vector<double> rEA_L(_n_elems);

  std::vector<double> rA_R(_n_elems);
  std::vector<double> ruA_R(_n_elems);
  std::vector<double> rEA_R(_n_elems);

  // transient
  unsigned int k = 1; // time step index
  while (_in_transient)
  {
    double max_wave_speed = 0;
    for (unsigned int i = 0; i < _n_elems; i++)
    {
      const double r = rA[0][i] / _A_elem[i];
      const double u = ruA[0][i] / rA[0][i];
      const double E = rEA[0][i] / rA[0][i];
      const double e = E - 0.5 * u * u;
      const double c = _eos.c_from_r_e(r, e);
      max_wave_speed = std::max(max_wave_speed, std::abs(u) + c);
    }
    const double dt_nominal = _time_step_sizer.computeTimeStepSize(max_wave_speed, _dx_min);

    // Check dt does not over-extend end time and flag end of transient
    const double dt = getTimeStepSizeAndUpdateTransientFlag(dt_nominal, t);

    const double CFL = Numerics::computeCFL(dt, max_wave_speed, _dx_min);

    std::cout << "Time step " << k << ": t = " << t + dt << " s, dt = " << dt
      << " s, CFL = " << CFL << std::endl;

    for (unsigned int s = 1; s < _n_stages+1; s++)
    {
      // compute solution on interfaces
      reconstructor.reconstructSolution(rA[s-1], ruA[s-1], rEA[s-1], _A_elem, _A_node, _x_elem, _x_node,
        rA_L, ruA_L, rEA_L, rA_R, ruA_R, rEA_R);

      // compute fluxes
      std::vector<std::vector<double>> f(_n_nodes, std::vector<double>(3, 0.0));
      {
        const unsigned int in = 0;
        const unsigned int ib = 0;
        const std::vector<double> U_b = {rA_L[ib], ruA_L[ib], rEA_L[ib]};
        f[in] = bc_left.computeFlux(U_b, _A_node[in]);
      }
      for (unsigned int in = 1; in < _n_nodes - 1; in++)
      {
        const unsigned int iL = in - 1;
        const unsigned int iR = in;
        const std::vector<double> U_L = {rA_R[iL], ruA_R[iL], rEA_R[iL]};
        const std::vector<double> U_R = {rA_L[iR], ruA_L[iR], rEA_L[iR]};
        f[in] = flux.computeFlux(U_L, U_R, _A_node[in]);
      }
      {
        const unsigned int in = _n_nodes - 1;
        const unsigned int ib = _n_elems - 1;
        const std::vector<double> U_b = {rA_R[ib], ruA_R[ib], rEA_R[ib]};
        f[in] = bc_right.computeFlux(U_b, _A_node[in]);
      }

      // solve
      for (unsigned int i = 0; i < _n_elems; i++)
      {
        const double r = rA[s-1][i] / _A_elem[i];
        const double u = ruA[s-1][i] / rA[s-1][i];
        const double E = rEA[s-1][i] / rA[s-1][i];
        const double e = E - 0.5 * u * u;
        const double p = _eos.p_from_r_e(r, e);

        rA[s][i] = _rk_b[s-1] * dt * (f[i][0] - f[i + 1][0]) / _dx;
        ruA[s][i] = _rk_b[s-1] * dt * (f[i][1] - f[i + 1][1] + p * (_A_node[i + 1] - _A_node[i])) / _dx;
        rEA[s][i] = _rk_b[s-1] * dt * (f[i][2] - f[i + 1][2]) / _dx;

        for (unsigned int k = 0; k <= s - 1; k++)
        {
          rA[s][i] += _rk_a[s-1][k] * rA[k][i];
          ruA[s][i] += _rk_a[s-1][k] * ruA[k][i];
          rEA[s][i] += _rk_a[s-1][k] * rEA[k][i];
        }
      }
    }

    rA[0] = rA[_n_stages];
    ruA[0] = ruA[_n_stages];
    rEA[0] = rEA[_n_stages];

    t += dt;
    k += 1;
  }

  // output

  std::vector<double> r(_n_elems);
  std::vector<double> u(_n_elems);
  std::vector<double> p(_n_elems);
  std::vector<double> T(_n_elems);

  for (unsigned int i = 0; i < _n_elems; i++)
  {
    r[i] = rA[0][i] / _A_elem[i];
    u[i] = ruA[0][i] / rA[0][i];
    const double E = rEA[0][i] / rA[0][i];
    const double e = E - 0.5 * u[i] * u[i];
    p[i] = _eos.p_from_r_e(r[i], e);
    T[i] = _eos.T_from_r_e(r[i], e);
  }

  CSVOutput csv_output("output.csv");
  csv_output.addOutput(_x_elem, "x");
  csv_output.addOutput(_A_elem, "A");
  csv_output.addOutput(r, "r");
  csv_output.addOutput(u, "u");
  csv_output.addOutput(p, "p");
  csv_output.addOutput(T, "T");
  csv_output.save();
}

std::vector<double> ExecutionerEuler1Phase::computeAreaNode() const
{
  std::vector<double> A_node(_n_nodes);
  for (unsigned int i = 0; i < _n_nodes; i++)
    A_node[i] = _A_fn.value(_x_node[i], 0.0);
  return A_node;
}

std::vector<double> ExecutionerEuler1Phase::computeAreaElem() const
{
  std::vector<double> A_elem(_n_elems);
  for (unsigned int i = 0; i < _n_elems; i++)
    A_elem[i] = 0.5 * (_A_node[i] + _A_node[i + 1]);
  return A_elem;
}
