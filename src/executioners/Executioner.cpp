#include "Executioner.h"
#include "Problem.h"
#include "RunParameters.h"
#include "utils.h"

Executioner::Executioner(
  const Problem & problem,
  const RunParameters & run_params)
  : _problem_base(problem),
    _run_params_base(run_params),

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
