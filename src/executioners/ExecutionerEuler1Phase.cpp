#include "ExecutionerEuler1Phase.h"
#include "ProblemEuler1Phase.h"
#include "RunParametersEuler1Phase.h"
#include "DoFHandlerEuler1Phase.h"
#include "BCEuler1Phase.h"
#include "ReconstructorEuler1Phase.h"
#include "FluxEuler1Phase.h"
#include "Function.h"
#include "EOS1Phase.h"
#include "CSVOutput.h"

ExecutionerEuler1Phase::ExecutionerEuler1Phase(
  const ProblemEuler1Phase & problem,
  const RunParametersEuler1Phase & run_params)
  : Executioner(problem, run_params, run_params.getNumberOfElements() * 3, 3),
    _problem(problem),
    _run_params(run_params),

    _dof_handler(run_params.getDoFHandler()),
    _eos(problem.getEOS()),
    _A_fn(problem.getAreaFunction()),
    _r_ic_fn(problem.getICDensity()),
    _u_ic_fn(problem.getICVelocity()),
    _p_ic_fn(problem.getICPressure()),
    _bc_left(_problem.getBCLeft()),
    _bc_right(_problem.getBCRight()),
    _flux(_run_params.getFlux()),
    _reconstructor(_run_params.getReconstructor()),

    _A_node(computeAreaNode()),
    _A_elem(computeAreaElem())
{
}

void ExecutionerEuler1Phase::initializeSolution(std::vector<double> & U) const
{
  for (unsigned int i = 0; i < _n_elems; i++)
  {
    const double r = _r_ic_fn.value(_x_elem[i], 0);
    const double u = _u_ic_fn.value(_x_elem[i], 0);
    const double p = _p_ic_fn.value(_x_elem[i], 0);

    const double e = _eos.e_from_p_r(p, r);
    const double E = e + 0.5 * u * u;

    const unsigned int i_rA = _dof_handler.elemIndex_rA(i);
    const unsigned int i_ruA = _dof_handler.elemIndex_ruA(i);
    const unsigned int i_rEA = _dof_handler.elemIndex_rEA(i);

    U[i_rA] = r * _A_elem[i];
    U[i_ruA] = r * u * _A_elem[i];
    U[i_rEA] = r * E * _A_elem[i];
  }
}

double ExecutionerEuler1Phase::computeMaxWaveSpeed(const std::vector<double> & U) const
{
  double max_wave_speed = 0;
  for (unsigned int i = 0; i < _n_elems; i++)
  {
    const unsigned int i_rA = _dof_handler.elemIndex_rA(i);
    const unsigned int i_ruA = _dof_handler.elemIndex_ruA(i);
    const unsigned int i_rEA = _dof_handler.elemIndex_rEA(i);

    const double r = U[i_rA] / _A_elem[i];
    const double u = U[i_ruA] / U[i_rA];
    const double E = U[i_rEA] / U[i_rA];
    const double e = E - 0.5 * u * u;
    const double c = _eos.c_from_r_e(r, e);
    max_wave_speed = std::max(max_wave_speed, std::abs(u) + c);
  }
  return max_wave_speed;
}

void ExecutionerEuler1Phase::computeSteadyStateResidual(
  const std::vector<double> & U, std::vector<double> & ss_rhs) const
{
  // compute solution on interfaces
  std::vector<double> U_L(_n_dofs);
  std::vector<double> U_R(_n_dofs);
  _reconstructor.reconstructSolution(U, _A_elem, _A_node, _x_elem, _x_node, U_L, U_R);

  // compute fluxes
  std::vector<std::vector<double>> f(_n_nodes, std::vector<double>(_n_vars, 0.0));
  {
    const unsigned int in = 0;
    const unsigned int ib = 0;
    const std::vector<double> U_b = _dof_handler.getElemSolutionVector(U_L, ib);
    f[in] = _bc_left.computeFlux(U_b, _A_node[in]);
  }
  for (unsigned int in = 1; in < _n_nodes - 1; in++)
  {
    const unsigned int iL = in - 1;
    const unsigned int iR = in;
    const std::vector<double> U_iL = _dof_handler.getElemSolutionVector(U_R, iL);
    const std::vector<double> U_iR = _dof_handler.getElemSolutionVector(U_L, iR);
    f[in] = _flux.computeFlux(U_iL, U_iR, _A_node[in]);
  }
  {
    const unsigned int in = _n_nodes - 1;
    const unsigned int ib = _n_elems - 1;
    const std::vector<double> U_b = _dof_handler.getElemSolutionVector(U_R, ib);
    f[in] = _bc_right.computeFlux(U_b, _A_node[in]);
  }

  for (unsigned int i = 0; i < _n_elems; i++)
  {
    const unsigned int i_rA = _dof_handler.elemIndex_rA(i);
    const unsigned int i_ruA = _dof_handler.elemIndex_ruA(i);
    const unsigned int i_rEA = _dof_handler.elemIndex_rEA(i);

    const double r = U[i_rA] / _A_elem[i];
    const double u = U[i_ruA] / U[i_rA];
    const double E = U[i_rEA] / U[i_rA];
    const double e = E - 0.5 * u * u;
    const double p = _eos.p_from_r_e(r, e);

    ss_rhs[i_rA] = (f[i][0] - f[i + 1][0]) / _dx;
    ss_rhs[i_ruA] = (f[i][1] - f[i + 1][1] + p * (_A_node[i + 1] - _A_node[i])) / _dx;
    ss_rhs[i_rEA] = (f[i][2] - f[i + 1][2]) / _dx;
  }
}

void ExecutionerEuler1Phase::outputSolution(const std::vector<double> & U) const
{
  std::vector<double> r(_n_elems);
  std::vector<double> u(_n_elems);
  std::vector<double> p(_n_elems);
  std::vector<double> T(_n_elems);

  for (unsigned int i = 0; i < _n_elems; i++)
  {
    const unsigned int i_rA = _dof_handler.elemIndex_rA(i);
    const unsigned int i_ruA = _dof_handler.elemIndex_ruA(i);
    const unsigned int i_rEA = _dof_handler.elemIndex_rEA(i);

    r[i] = U[i_rA] / _A_elem[i];
    u[i] = U[i_ruA] / U[i_rA];
    const double E = U[i_rEA] / U[i_rA];
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
