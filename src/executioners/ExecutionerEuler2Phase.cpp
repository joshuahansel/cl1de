#include "ExecutionerEuler2Phase.h"
#include "ProblemEuler2Phase.h"
#include "RunParametersEuler2Phase.h"
#include "DoFHandlerEuler2Phase.h"
#include "ICsEuler2Phase.h"
#include "BCEuler2Phase.h"
#include "ReconstructorEuler2Phase.h"
#include "FluxEuler2Phase.h"
#include "Function.h"
#include "EOS1Phase.h"
#include "EOS1PhaseStiffenedGas.h"
#include "CSVOutput.h"

#include <cmath>
#include <iostream>

ExecutionerEuler2Phase::ExecutionerEuler2Phase(
  const ProblemEuler2Phase & problem,
  const RunParametersEuler2Phase & run_params)
  : Executioner(problem, run_params, run_params.getNumberOfElements() * 7, 7),
    _problem(problem),
    _run_params(run_params),

    _dof_handler(run_params.getDoFHandler()),
    _eos_liq(problem.getEOSLiquid()),
    _eos_vap(problem.getEOSVapor()),
    _A_fn(problem.getAreaFunction()),
    _ics(problem.getICs()),
    _bc_left(_problem.getBCLeft()),
    _bc_right(_problem.getBCRight()),
    _gravity(_problem.getGravity()),
    _a_int(_problem.getInterfacialAreaDensity()),
    _flux(_run_params.getFlux()),
    _reconstructor(_run_params.getReconstructor()),

    _use_finite_velocity_relaxation(_run_params.getUseFiniteVelocityRelaxation()),
    _use_finite_pressure_relaxation(_run_params.getUseFinitePressureRelaxation()),

    _use_infinite_velocity_relaxation(_run_params.getUseInfiniteVelocityRelaxation()),
    _use_infinite_pressure_relaxation(_run_params.getUseInfinitePressureRelaxation()),

    _n_source_substeps(_run_params.getNumberOfSourceSubsteps()),

    _A_node(computeAreaNode()),
    _A_elem(computeAreaElem())
{
}

void ExecutionerEuler2Phase::initializeSolution(std::vector<double> & U) const
{
  for (unsigned int i = 0; i < _n_elems; i++)
  {
    const unsigned int i_aA_liq = _dof_handler.elemIndex_aA_liq(i);
    const unsigned int i_arA_liq = _dof_handler.elemIndex_arA_liq(i);
    const unsigned int i_aruA_liq = _dof_handler.elemIndex_aruA_liq(i);
    const unsigned int i_arEA_liq = _dof_handler.elemIndex_arEA_liq(i);
    const unsigned int i_arA_vap = _dof_handler.elemIndex_arA_vap(i);
    const unsigned int i_aruA_vap = _dof_handler.elemIndex_aruA_vap(i);
    const unsigned int i_arEA_vap = _dof_handler.elemIndex_arEA_vap(i);

    _ics.computeICs(
      _x_elem[i],
      _A_elem[i],
      U[i_aA_liq],
      U[i_arA_liq],
      U[i_aruA_liq],
      U[i_arEA_liq],
      U[i_arA_vap],
      U[i_aruA_vap],
      U[i_arEA_vap]);
  }
}

double ExecutionerEuler2Phase::computeMaxWaveSpeed(const std::vector<double> & U) const
{
  double max_wave_speed = 0;
  for (unsigned int i = 0; i < _n_elems; i++)
  {
    const unsigned int i_aA_liq = _dof_handler.elemIndex_aA_liq(i);
    const unsigned int i_arA_liq = _dof_handler.elemIndex_arA_liq(i);
    const unsigned int i_aruA_liq = _dof_handler.elemIndex_aruA_liq(i);
    const unsigned int i_arEA_liq = _dof_handler.elemIndex_arEA_liq(i);
    const unsigned int i_arA_vap = _dof_handler.elemIndex_arA_vap(i);
    const unsigned int i_aruA_vap = _dof_handler.elemIndex_aruA_vap(i);
    const unsigned int i_arEA_vap = _dof_handler.elemIndex_arEA_vap(i);

    const double A = _A_elem[i];

    const double a_liq = U[i_aA_liq] / A;
    const double r_liq = U[i_arA_liq] / U[i_aA_liq];
    const double u_liq = U[i_aruA_liq] / U[i_arA_liq];
    const double E_liq = U[i_arEA_liq] / U[i_arA_liq];

    const double a_vap = 1.0 - a_liq;
    const double r_vap = U[i_arA_vap] / (a_vap * A);
    const double u_vap = U[i_aruA_vap] / U[i_arA_vap];
    const double E_vap = U[i_arEA_vap] / U[i_arA_vap];

    const double e_liq = E_liq - 0.5 * u_liq * u_liq;
    const double c_liq = _eos_liq.c_from_r_e(r_liq, e_liq);
    const double e_vap = E_vap - 0.5 * u_vap * u_vap;
    const double c_vap = _eos_vap.c_from_r_e(r_vap, e_vap);

    const double wave_speed = std::max(std::abs(u_liq) + c_liq, std::abs(u_vap) + c_vap);
    max_wave_speed = std::max(max_wave_speed, wave_speed);
  }
  return max_wave_speed;
}

void ExecutionerEuler2Phase::computeSteadyStateResidual(
  const std::vector<double> & U, std::vector<double> & ss_rhs) const
{
  // compute solution on interfaces
  std::vector<double> U_L(_n_dofs);
  std::vector<double> U_R(_n_dofs);
  _reconstructor.reconstructSolution(U, _A_elem, _A_node, _x_elem, _x_node, U_L, U_R);

  // compute fluxes
  std::vector<std::vector<double>> f_L(_n_elems, std::vector<double>(_n_vars, 0.0));
  std::vector<std::vector<double>> f_R(_n_elems, std::vector<double>(_n_vars, 0.0));
  {
    const unsigned int in = 0;
    const unsigned int ib = 0;
    const std::vector<double> U_b = _dof_handler.getElemSolutionVector(U_L, ib);
    f_L[ib] = _bc_left.computeFlux(U_b, _A_node[in]);
  }
  for (unsigned int in = 1; in < _n_nodes - 1; in++)
  {
    const unsigned int iL = in - 1;
    const unsigned int iR = in;
    const std::vector<double> U_iL = _dof_handler.getElemSolutionVector(U_R, iL);
    const std::vector<double> U_iR = _dof_handler.getElemSolutionVector(U_L, iR);

    _flux.computeFlux(U_iL, U_iR, _A_node[in], f_R[iL], f_L[iR]);
  }
  {
    const unsigned int in = _n_nodes - 1;
    const unsigned int ib = _n_elems - 1;
    const std::vector<double> U_b = _dof_handler.getElemSolutionVector(U_R, ib);
    f_R[ib] = _bc_right.computeFlux(U_b, _A_node[in]);
  }

  for (unsigned int i = 0; i < _n_elems; i++)
  {
    const unsigned int i_aA_liq = _dof_handler.elemIndex_aA_liq(i);
    const unsigned int i_arA_liq = _dof_handler.elemIndex_arA_liq(i);
    const unsigned int i_aruA_liq = _dof_handler.elemIndex_aruA_liq(i);
    const unsigned int i_arEA_liq = _dof_handler.elemIndex_arEA_liq(i);
    const unsigned int i_arA_vap = _dof_handler.elemIndex_arA_vap(i);
    const unsigned int i_aruA_vap = _dof_handler.elemIndex_aruA_vap(i);
    const unsigned int i_arEA_vap = _dof_handler.elemIndex_arEA_vap(i);

    const double A = _A_elem[i];

    const double a_liq = U[i_aA_liq] / A;
    const double r_liq = U[i_arA_liq] / U[i_aA_liq];
    const double u_liq = U[i_aruA_liq] / U[i_arA_liq];
    const double E_liq = U[i_arEA_liq] / U[i_arA_liq];

    const double a_vap = 1.0 - a_liq;
    const double r_vap = U[i_arA_vap] / (a_vap * A);
    const double u_vap = U[i_aruA_vap] / U[i_arA_vap];
    const double E_vap = U[i_arEA_vap] / U[i_arA_vap];

    const double e_liq = E_liq - 0.5 * u_liq * u_liq;
    const double e_vap = E_vap - 0.5 * u_vap * u_vap;

    const double p_liq = _eos_liq.p_from_r_e(r_liq, e_liq);
    const double p_vap = _eos_vap.p_from_r_e(r_vap, e_vap);

    ss_rhs[i_aA_liq]   = (f_L[i][0] - f_R[i][0]) / _dx;
    ss_rhs[i_arA_liq]  = (f_L[i][1] - f_R[i][1]) / _dx;
    ss_rhs[i_aruA_liq] = (f_L[i][2] - f_R[i][2]
      + a_liq * p_liq * (_A_node[i + 1] - _A_node[i])) / _dx
      + a_liq * r_liq * _gravity * _A_elem[i];
    ss_rhs[i_arEA_liq] = (f_L[i][3] - f_R[i][3]) / _dx
      + a_liq * r_liq * u_liq * _gravity * _A_elem[i];
    ss_rhs[i_arA_vap]  = (f_L[i][4] - f_R[i][4]) / _dx;
    ss_rhs[i_aruA_vap] = (f_L[i][5] - f_R[i][5]
      + a_vap * p_vap * (_A_node[i + 1] - _A_node[i])) / _dx
      + a_vap * r_vap * _gravity * _A_elem[i];
    ss_rhs[i_arEA_vap] = (f_L[i][6] - f_R[i][6]) / _dx
      + a_vap * r_vap * u_vap * _gravity * _A_elem[i];
  }
}

void ExecutionerEuler2Phase::performPostStep(std::vector<double> & U, const double & dt) const
{
  double relax_u, relax_p, u_int_bar, p_int_bar;
  if (_use_finite_velocity_relaxation || _use_finite_pressure_relaxation)
  {
    const double sub_dt = dt / _n_source_substeps;
    for (unsigned int n = 0; n < _n_source_substeps; n++)
      for (unsigned int i = 0; i < _n_elems; i++)
      {
        const unsigned int i_aA_liq = _dof_handler.elemIndex_aA_liq(i);
        const unsigned int i_arA_liq = _dof_handler.elemIndex_arA_liq(i);
        const unsigned int i_aruA_liq = _dof_handler.elemIndex_aruA_liq(i);
        const unsigned int i_arEA_liq = _dof_handler.elemIndex_arEA_liq(i);
        const unsigned int i_arA_vap = _dof_handler.elemIndex_arA_vap(i);
        const unsigned int i_aruA_vap = _dof_handler.elemIndex_aruA_vap(i);
        const unsigned int i_arEA_vap = _dof_handler.elemIndex_arEA_vap(i);

        const double A = _A_elem[i];

        const double a_liq = U[i_aA_liq] / A;
        const double r_liq = U[i_arA_liq] / U[i_aA_liq];
        const double u_liq = U[i_aruA_liq] / U[i_arA_liq];
        const double E_liq = U[i_arEA_liq] / U[i_arA_liq];

        const double a_vap = 1.0 - a_liq;
        const double r_vap = U[i_arA_vap] / (a_vap * A);
        const double u_vap = U[i_aruA_vap] / U[i_arA_vap];
        const double E_vap = U[i_arEA_vap] / U[i_arA_vap];

        const double e_liq = E_liq - 0.5 * u_liq * u_liq;
        const double e_vap = E_vap - 0.5 * u_vap * u_vap;

        const double p_liq = _eos_liq.p_from_r_e(r_liq, e_liq);
        const double p_vap = _eos_vap.p_from_r_e(r_vap, e_vap);

        const double c_liq = _eos_liq.c_from_r_e(r_liq, e_liq);
        const double c_vap = _eos_vap.c_from_r_e(r_vap, e_vap);

        const double Z_liq = r_liq * c_liq;
        const double Z_vap = r_vap * c_vap;

        u_int_bar = (Z_liq * u_liq + Z_vap * u_vap) / (Z_liq + Z_vap);
        p_int_bar = (Z_liq * p_vap + Z_vap * p_liq) / (Z_liq + Z_vap);

        relax_p = _a_int / (Z_liq + Z_vap);
        relax_u = 0.5 * relax_p * Z_liq * Z_vap;

        double dU_aA_liq = 0, dU_aruA_liq = 0, dU_arEA_liq = 0, dU_aruA_vap = 0, dU_arEA_vap = 0;
        if (_use_finite_velocity_relaxation)
        {
          dU_aruA_liq += sub_dt * relax_u * (u_vap - u_liq) * _A_elem[i];
          dU_arEA_liq += sub_dt * relax_u * u_int_bar * (u_vap - u_liq) * _A_elem[i];
          dU_aruA_vap -= sub_dt * relax_u * (u_vap - u_liq) * _A_elem[i];
          dU_arEA_vap -= sub_dt * relax_u * u_int_bar * (u_vap - u_liq) * _A_elem[i];
        }
        if (_use_finite_pressure_relaxation)
        {
          dU_aA_liq   += sub_dt * relax_p * (p_liq - p_vap) * _A_elem[i];
          dU_arEA_liq -= sub_dt * relax_p * p_int_bar * (p_liq - p_vap) * _A_elem[i];
          dU_arEA_vap += sub_dt * relax_p * p_int_bar * (p_liq - p_vap) * _A_elem[i];
        }

        U[i_aA_liq] += dU_aA_liq;
        U[i_aruA_liq] += dU_aruA_liq;
        U[i_arEA_liq] += dU_arEA_liq;
        U[i_aruA_vap] += dU_aruA_vap;
        U[i_arEA_vap] += dU_arEA_vap;
      }
  }

  if (_use_infinite_velocity_relaxation)
  {
    for (unsigned int i = 0; i < _n_elems; i++)
    {
      const unsigned int i_aA_liq = _dof_handler.elemIndex_aA_liq(i);
      const unsigned int i_arA_liq = _dof_handler.elemIndex_arA_liq(i);
      const unsigned int i_aruA_liq = _dof_handler.elemIndex_aruA_liq(i);
      const unsigned int i_arEA_liq = _dof_handler.elemIndex_arEA_liq(i);
      const unsigned int i_arA_vap = _dof_handler.elemIndex_arA_vap(i);
      const unsigned int i_aruA_vap = _dof_handler.elemIndex_aruA_vap(i);
      const unsigned int i_arEA_vap = _dof_handler.elemIndex_arEA_vap(i);

      const double u_liq = U[i_aruA_liq] / U[i_arA_liq];
      const double E_liq = U[i_arEA_liq] / U[i_arA_liq];
      const double u_vap = U[i_aruA_vap] / U[i_arA_vap];
      const double E_vap = U[i_arEA_vap] / U[i_arA_vap];

      const double x_liq = U[i_arA_liq] / (U[i_arA_liq] + U[i_arA_vap]);
      const double x_vap = U[i_arA_vap] / (U[i_arA_liq] + U[i_arA_vap]);
      const double u_rel = x_liq + u_liq + x_vap + u_vap;
      const double E_liq_rel = E_liq + u_rel * (u_rel - u_liq);
      const double E_vap_rel = E_vap + u_rel * (u_rel - u_vap);

      U[i_aruA_liq] = U[i_arA_liq] * u_rel;
      U[i_aruA_vap] = U[i_arA_vap] * u_rel;
      U[i_arEA_liq] = U[i_arA_liq] * E_liq_rel;
      U[i_arEA_vap] = U[i_arA_vap] * E_vap_rel;
    }
  }

  if (_use_infinite_pressure_relaxation)
  {
    const auto & eos_liq = dynamic_cast<const EOS1PhaseStiffenedGas &>(_eos_liq);
    const auto & eos_vap = dynamic_cast<const EOS1PhaseStiffenedGas &>(_eos_vap);

    const double & gamma_liq = eos_liq.getGamma();
    const double & gamma_vap = eos_vap.getGamma();
    const double & p_inf_liq = eos_liq.getReferencePressure();
    const double & p_inf_vap = eos_vap.getReferencePressure();

    for (unsigned int i = 0; i < _n_elems; i++)
    {
      const unsigned int i_aA_liq = _dof_handler.elemIndex_aA_liq(i);
      const unsigned int i_arA_liq = _dof_handler.elemIndex_arA_liq(i);
      const unsigned int i_aruA_liq = _dof_handler.elemIndex_aruA_liq(i);
      const unsigned int i_arEA_liq = _dof_handler.elemIndex_arEA_liq(i);
      const unsigned int i_arA_vap = _dof_handler.elemIndex_arA_vap(i);
      const unsigned int i_aruA_vap = _dof_handler.elemIndex_aruA_vap(i);
      const unsigned int i_arEA_vap = _dof_handler.elemIndex_arEA_vap(i);

      const double A = _A_elem[i];

      const double a_liq = U[i_aA_liq] / A;
      const double r_liq = U[i_arA_liq] / U[i_aA_liq];
      const double u_liq = U[i_aruA_liq] / U[i_arA_liq];
      const double E_liq = U[i_arEA_liq] / U[i_arA_liq];

      const double a_vap = 1.0 - a_liq;
      const double r_vap = U[i_arA_vap] / (a_vap * A);
      const double u_vap = U[i_aruA_vap] / U[i_arA_vap];
      const double E_vap = U[i_arEA_vap] / U[i_arA_vap];

      const double v_liq = 1.0 / r_liq;
      const double v_vap = 1.0 / r_vap;

      const double e_liq = E_liq - 0.5 * u_liq * u_liq;
      const double e_vap = E_vap - 0.5 * u_vap * u_vap;

      const double p_liq = _eos_liq.p_from_r_e(r_liq, e_liq);
      const double p_vap = _eos_vap.p_from_r_e(r_vap, e_vap);

      const double A_liq = a_liq / gamma_liq * (p_liq + p_inf_liq)
        / (a_liq / gamma_liq + a_vap / gamma_vap);
      const double A_vap = a_vap / gamma_vap * (p_vap + p_inf_vap)
        / (a_liq / gamma_liq + a_vap / gamma_vap);
      // const double p_rel = 0.5 * (A_liq + A_vap - p_inf_liq - p_inf_vap)
      //   + std::sqrt(0.25 * std::pow(A_vap - A_liq - p_inf_vap + p_inf_liq, 2) + A_liq * A_vap);
      const double p_rel = 0.5 * (A_liq + A_vap - p_inf_liq - p_inf_vap)
        + std::sqrt(0.25 * std::pow(-A_vap - A_liq + p_inf_vap + p_inf_liq, 2)
        - p_inf_liq * p_inf_vap + A_liq * p_inf_vap + A_vap * p_inf_liq);
      const double a_liq_rel = a_liq / gamma_liq * (gamma_liq - 1.0 + (p_liq + p_inf_liq) / (p_rel + p_inf_liq));
      const double a_vap_rel = 1.0 - a_liq_rel;
      const double r_liq_rel = a_liq / a_liq_rel * r_liq;
      const double r_vap_rel = a_vap / a_vap_rel * r_vap;
      const double v_liq_rel = 1.0 / r_liq_rel;
      const double v_vap_rel = 1.0 / r_vap_rel;
      const double E_liq_rel = E_liq - p_rel * (v_liq_rel - v_liq);
      const double E_vap_rel = E_vap - p_rel * (v_vap_rel - v_vap);

      const double arA_liq_rel = a_liq_rel * r_liq_rel * A;
      const double arA_vap_rel = a_vap_rel * r_vap_rel * A;

      U[i_aA_liq]   = a_liq_rel * A;
      U[i_arA_liq]  = arA_liq_rel;
      U[i_aruA_liq] = arA_liq_rel * u_liq;
      U[i_arEA_liq] = arA_liq_rel * E_liq_rel;
      U[i_arA_vap]  = arA_vap_rel;
      U[i_aruA_vap] = arA_vap_rel * u_vap;
      U[i_arEA_vap] = arA_vap_rel * E_vap_rel;
    }
  }
}

void ExecutionerEuler2Phase::outputSolution(const std::vector<double> & U) const
{
  std::vector<double> a_liq(_n_elems);
  std::vector<double> r_liq(_n_elems);
  std::vector<double> u_liq(_n_elems);
  std::vector<double> p_liq(_n_elems);
  std::vector<double> T_liq(_n_elems);
  std::vector<double> a_vap(_n_elems);
  std::vector<double> r_vap(_n_elems);
  std::vector<double> u_vap(_n_elems);
  std::vector<double> p_vap(_n_elems);
  std::vector<double> T_vap(_n_elems);

  std::vector<double> r_mix(_n_elems);
  std::vector<double> u_mix(_n_elems);
  std::vector<double> p_mix(_n_elems);

  for (unsigned int i = 0; i < _n_elems; i++)
  {
    const unsigned int i_aA_liq = _dof_handler.elemIndex_aA_liq(i);
    const unsigned int i_arA_liq = _dof_handler.elemIndex_arA_liq(i);
    const unsigned int i_aruA_liq = _dof_handler.elemIndex_aruA_liq(i);
    const unsigned int i_arEA_liq = _dof_handler.elemIndex_arEA_liq(i);
    const unsigned int i_arA_vap = _dof_handler.elemIndex_arA_vap(i);
    const unsigned int i_aruA_vap = _dof_handler.elemIndex_aruA_vap(i);
    const unsigned int i_arEA_vap = _dof_handler.elemIndex_arEA_vap(i);

    const double A = _A_elem[i];

    a_liq[i] = U[i_aA_liq] / A;
    r_liq[i] = U[i_arA_liq] / U[i_aA_liq];
    u_liq[i] = U[i_aruA_liq] / U[i_arA_liq];
    const double E_liq = U[i_arEA_liq] / U[i_arA_liq];
    const double e_liq = E_liq - 0.5 * u_liq[i] * u_liq[i];
    p_liq[i] = _eos_liq.p_from_r_e(r_liq[i], e_liq);
    T_liq[i] = _eos_liq.T_from_r_e(r_liq[i], e_liq);

    a_vap[i] = 1.0 - a_liq[i];
    r_vap[i] = U[i_arA_vap] / (a_vap[i] * A);
    u_vap[i] = U[i_aruA_vap] / U[i_arA_vap];
    const double E_vap = U[i_arEA_vap] / U[i_arA_vap];
    const double e_vap = E_vap - 0.5 * u_vap[i] * u_vap[i];
    p_vap[i] = _eos_vap.p_from_r_e(r_vap[i], e_vap);
    T_vap[i] = _eos_vap.T_from_r_e(r_vap[i], e_vap);

    r_mix[i] = a_liq[i] * r_liq[i] + a_vap[i] * r_vap[i];
    u_mix[i] = a_liq[i] * u_liq[i] + a_vap[i] * u_vap[i];
    p_mix[i] = a_liq[i] * p_liq[i] + a_vap[i] * p_vap[i];
  }

  CSVOutput csv_output("output.csv");
  csv_output.addOutput(_x_elem, "x");
  csv_output.addOutput(_A_elem, "A");
  csv_output.addOutput(a_liq, "a_liq");
  csv_output.addOutput(r_liq, "r_liq");
  csv_output.addOutput(u_liq, "u_liq");
  csv_output.addOutput(p_liq, "p_liq");
  csv_output.addOutput(T_liq, "T_liq");
  csv_output.addOutput(a_vap, "a_vap");
  csv_output.addOutput(r_vap, "r_vap");
  csv_output.addOutput(u_vap, "u_vap");
  csv_output.addOutput(p_vap, "p_vap");
  csv_output.addOutput(T_vap, "T_vap");
  csv_output.addOutput(r_mix, "r_mix");
  csv_output.addOutput(u_mix, "u_mix");
  csv_output.addOutput(p_mix, "p_mix");
  csv_output.save();
}

std::vector<double> ExecutionerEuler2Phase::computeAreaNode() const
{
  std::vector<double> A_node(_n_nodes);
  for (unsigned int i = 0; i < _n_nodes; i++)
    A_node[i] = _A_fn.value(_x_node[i], 0.0);
  return A_node;
}

std::vector<double> ExecutionerEuler2Phase::computeAreaElem() const
{
  std::vector<double> A_elem(_n_elems);
  for (unsigned int i = 0; i < _n_elems; i++)
    A_elem[i] = 0.5 * (_A_node[i] + _A_node[i + 1]);
  return A_elem;
}
