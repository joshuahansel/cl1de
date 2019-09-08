#include "Executioner.h"
#include "Problem.h"
#include "RunParameters.h"
#include "TimeStepSizer.h"
#include "utils.h"
#include "UtilsNumerics.h"

#include <iostream>
#include <vector>

Executioner::Executioner(
  const Problem & problem,
  const RunParameters & run_params,
  const unsigned int & n_dofs,
  const unsigned int & n_vars)
  : _problem_base(problem),
    _run_params_base(run_params),
    _n_dofs(n_dofs),
    _n_vars(n_vars),

    _n_elems(_run_params_base.getNumberOfElements()),
    _n_nodes(_n_elems + 1),
    _dx(computeElemSize()),
    _dx_min(_dx),
    _x_elem(computeElemPositions()),
    _x_node(computeNodePositions()),

    _t_end(getEndTime()),
    _n_stages(_run_params_base.getNumberOfStages()),
    _rk_a(initializeRungeKuttaSolutionCoefs()),
    _rk_b(initializeRungeKuttaTimeStepCoefs()),

    _time_step_sizer(run_params.getTimeStepSizer()),

    _in_transient(true)
{
}

void Executioner::run()
{
  std::vector<std::vector<double>> U(_n_stages + 1, std::vector<double>(_n_dofs, 0));

  // initialize
  initializeSolution(U[0]);

  std::vector<double> ss_rhs(_n_dofs, 0.0);

  // transient
  double t = 0.0;
  unsigned int k = 1; // time step index
  while (_in_transient)
  {
    const double max_wave_speed = computeMaxWaveSpeed(U[0]);
    const double dt_nominal = _time_step_sizer.computeTimeStepSize(max_wave_speed, _dx_min);

    // Check dt does not over-extend end time and flag end of transient
    const double dt = getTimeStepSizeAndUpdateTransientFlag(dt_nominal, t);

    const double CFL = Numerics::computeCFL(dt, max_wave_speed, _dx_min);

    std::cout << "Time step " << k << ": t = " << t + dt << " s, dt = " << dt
      << " s, CFL = " << CFL << std::endl;

    for (unsigned int s = 1; s < _n_stages+1; s++)
    {
      computeSteadyStateResidual(U[s-1], ss_rhs);

      for (unsigned int i = 0; i < _n_dofs; i++)
      {
        U[s][i] = _rk_b[s-1] * dt * ss_rhs[i];
        for (unsigned int k = 0; k <= s - 1; k++)
          U[s][i] += _rk_a[s-1][k] * U[k][i];
      }
    }

    performPostStep(U[_n_stages], dt);

    U[0] = U[_n_stages];

    t += dt;
    k += 1;
  }

  // output
  outputSolution(U[0]);
}

void Executioner::performPostStep(std::vector<double> & /*U*/, const double & /*dt*/) const
{
}

double Executioner::computeElemSize() const
{
  const double & x_min = _problem_base.getXMin();
  const double & x_max = _problem_base.getXMax();
  return (x_max - x_min) / _n_elems;
}

std::vector<double> Executioner::computeElemPositions() const
{
  const double & x_min = _problem_base.getXMin();

  std::vector<double> x(_n_elems);
  for (unsigned int i = 0; i < _n_elems; i++)
    x[i] = x_min + (i + 0.5) * _dx;

  return x;
}

std::vector<double> Executioner::computeNodePositions() const
{
  const double & x_min = _problem_base.getXMin();

  std::vector<double> x(_n_nodes);
  x[0] = x_min;
  for (unsigned int i = 0; i < _n_elems; i++)
    x[i + 1] = x[i] + _dx;

  return x;
}

double Executioner::getEndTime() const
{
  if (_run_params_base.specifiedEndTime())
    return _run_params_base.getEndTime();
  else if (_problem_base.hasDefaultEndTime())
    return _problem_base.getDefaultEndTime();
  else
    throwError("No end time was specified, and the problem has no default end time.",
      __PRETTY_FUNCTION__);
}

double Executioner::getTimeStepSizeAndUpdateTransientFlag(const double & dt_nominal, const double & t)
{
  if (t + dt_nominal >= _t_end - 1e-15)
  {
    _in_transient = false;
    return _t_end - t;
  }
  else
    return dt_nominal;
}

std::vector<std::vector<double>> Executioner::initializeRungeKuttaSolutionCoefs() const
{
  std::vector<std::vector<double>> a(_n_stages, std::vector<double>(_n_stages, 0.0));
  if (_n_stages == 1)
  {
    a[0][0] = 1.0;
  }
  else if (_n_stages == 2)
  {
    a[0][0] = 1.0;
    a[1][0] = 0.5;
    a[1][1] = 0.5;
    // a[0][0] = 1.0;
    // a[1][0] = 1.0;
  }
  else if (_n_stages == 3)
  {
    a[0][0] = 1.0;
    a[1][0] = 0.75;
    a[1][1] = 0.25;
    a[2][0] = 1.0 / 3.0;
    a[2][1] = 0.0;
    a[2][2] = 2.0 / 3.0;
  }
  else
    throwError("Invalid number of time integrator stages.");

  return a;
}

std::vector<double> Executioner::initializeRungeKuttaTimeStepCoefs() const
{
  std::vector<double> b(_n_stages, 0.0);
  if (_n_stages == 1)
  {
    b[0] = 1.0;
  }
  else if (_n_stages == 2)
  {
    b[0] = 1.0;
    b[1] = 0.5;
    // b[0] = 0.5;
    // b[1] = 1.0;
  }
  else if (_n_stages == 3)
  {
    b[0] = 1.0;
    b[1] = 0.25;
    b[2] = 2.0 / 3.0;
  }
  else
    throwError("Invalid number of time integrator stages.");

  return b;
}
