#ifndef Executioner_H
#define Executioner_H

#include <vector>

class Problem;
class RunParameters;
class TimeStepSizer;

class Executioner
{
public:
  Executioner(
    const Problem & problem,
    const RunParameters & run_params,
    const unsigned int & n_dofs,
    const unsigned int & n_vars);

  void run();

protected:
  double computeElemSize() const;
  std::vector<double> computeElemPositions() const;
  std::vector<double> computeNodePositions() const;
  double getEndTime() const;
  double getTimeStepSizeAndUpdateTransientFlag(const double & dt_nominal, const double & t);
  std::vector<std::vector<double>> initializeRungeKuttaSolutionCoefs() const;
  std::vector<double> initializeRungeKuttaTimeStepCoefs() const;

  virtual void initializeSolution(std::vector<double> & U) const = 0;
  virtual double computeMaxWaveSpeed(const std::vector<double> & U) const = 0;
  virtual void computeSteadyStateResidual(
    const std::vector<double> & U, std::vector<double> & ss_rhs) const = 0;
  virtual void outputSolution(const std::vector<double> & U) const = 0;

  const Problem & _problem_base;
  const RunParameters & _run_params_base;
  const unsigned int _n_vars;
  const unsigned int _n_dofs;

  const unsigned int & _n_elems;
  const unsigned int _n_nodes;
  const double _dx;
  const double _dx_min;
  const std::vector<double> _x_elem;
  const std::vector<double> _x_node;

  const double _t_end;
  const unsigned int _n_stages;
  const std::vector<std::vector<double>> _rk_a;
  const std::vector<double> _rk_b;

  const TimeStepSizer & _time_step_sizer;

  bool _in_transient;
};

#endif /* Executioner_H */
