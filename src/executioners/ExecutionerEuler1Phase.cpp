#include "ExecutionerEuler1Phase.h"
#include "ProblemEuler1Phase.h"
#include "RunParametersEuler1Phase.h"
#include "DoFHandlerEuler1Phase.h"
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

    _dof_handler(run_params.getDoFHandler()),
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
  const unsigned int & n_dofs = _dof_handler.nDoFs();

  std::vector<std::vector<double>> U(_n_stages + 1, std::vector<double>(n_dofs, 0));

  // initialize
  double t = 0.0;

  for (unsigned int i = 0; i < _n_elems; i++)
  {
    const double r = _r_ic_fn.value(_x_elem[i], t);
    const double u = _u_ic_fn.value(_x_elem[i], t);
    const double p = _p_ic_fn.value(_x_elem[i], t);

    const double e = _eos.e_from_p_r(p, r);
    const double E = e + 0.5 * u * u;

    const unsigned int i_rA = _dof_handler.elemIndex_rA(i);
    const unsigned int i_ruA = _dof_handler.elemIndex_ruA(i);
    const unsigned int i_rEA = _dof_handler.elemIndex_rEA(i);

    U[0][i_rA] = r * _A_elem[i];
    U[0][i_ruA] = r * u * _A_elem[i];
    U[0][i_rEA] = r * E * _A_elem[i];
  }

  const BCEuler1Phase & bc_left = _problem.getBCLeft();
  const BCEuler1Phase & bc_right = _problem.getBCRight();

  const ReconstructorEuler1Phase & reconstructor = _run_params.getReconstructor();
  const FluxEuler1Phase & flux = _run_params.getFlux();

  std::vector<double> U_L(n_dofs);
  std::vector<double> U_R(n_dofs);

  const unsigned int & n_vars = _dof_handler.nVars();

  std::vector<double> ss_rhs(n_dofs, 0.0);

  // transient
  unsigned int k = 1; // time step index
  while (_in_transient)
  {
    double max_wave_speed = 0;
    for (unsigned int i = 0; i < _n_elems; i++)
    {
      const unsigned int i_rA = _dof_handler.elemIndex_rA(i);
      const unsigned int i_ruA = _dof_handler.elemIndex_ruA(i);
      const unsigned int i_rEA = _dof_handler.elemIndex_rEA(i);

      const double r = U[0][i_rA] / _A_elem[i];
      const double u = U[0][i_ruA] / U[0][i_rA];
      const double E = U[0][i_rEA] / U[0][i_rA];
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
      reconstructor.reconstructSolution(U[s-1], _A_elem, _A_node, _x_elem, _x_node, U_L, U_R);

      // compute fluxes
      std::vector<std::vector<double>> f(_n_nodes, std::vector<double>(n_vars, 0.0));
      {
        const unsigned int in = 0;
        const unsigned int ib = 0;
        const std::vector<double> U_b = _dof_handler.getElemSolutionVector(U_L, ib);
        f[in] = bc_left.computeFlux(U_b, _A_node[in]);
      }
      for (unsigned int in = 1; in < _n_nodes - 1; in++)
      {
        const unsigned int iL = in - 1;
        const unsigned int iR = in;
        const std::vector<double> U_iL = _dof_handler.getElemSolutionVector(U_R, iL);
        const std::vector<double> U_iR = _dof_handler.getElemSolutionVector(U_L, iR);
        f[in] = flux.computeFlux(U_iL, U_iR, _A_node[in]);
      }
      {
        const unsigned int in = _n_nodes - 1;
        const unsigned int ib = _n_elems - 1;
        const std::vector<double> U_b = _dof_handler.getElemSolutionVector(U_R, ib);
        f[in] = bc_right.computeFlux(U_b, _A_node[in]);
      }

      // solve
      for (unsigned int i = 0; i < _n_elems; i++)
      {
        const unsigned int i_rA = _dof_handler.elemIndex_rA(i);
        const unsigned int i_ruA = _dof_handler.elemIndex_ruA(i);
        const unsigned int i_rEA = _dof_handler.elemIndex_rEA(i);

        const double r = U[s-1][i_rA] / _A_elem[i];
        const double u = U[s-1][i_ruA] / U[s-1][i_rA];
        const double E = U[s-1][i_rEA] / U[s-1][i_rA];
        const double e = E - 0.5 * u * u;
        const double p = _eos.p_from_r_e(r, e);

        ss_rhs[i_rA] = (f[i][0] - f[i + 1][0]) / _dx;
        ss_rhs[i_ruA] = (f[i][1] - f[i + 1][1] + p * (_A_node[i + 1] - _A_node[i])) / _dx;
        ss_rhs[i_rEA] = (f[i][2] - f[i + 1][2]) / _dx;
      }

      for (unsigned int i = 0; i < n_dofs; i++)
      {
        U[s][i] = _rk_b[s-1] * dt * ss_rhs[i];
        for (unsigned int k = 0; k <= s - 1; k++)
          U[s][i] += _rk_a[s-1][k] * U[k][i];
      }
    }

    U[0] = U[_n_stages];

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
    const unsigned int i_rA = _dof_handler.elemIndex_rA(i);
    const unsigned int i_ruA = _dof_handler.elemIndex_ruA(i);
    const unsigned int i_rEA = _dof_handler.elemIndex_rEA(i);

    r[i] = U[0][i_rA] / _A_elem[i];
    u[i] = U[0][i_ruA] / U[0][i_rA];
    const double E = U[0][i_rEA] / U[0][i_rA];
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
