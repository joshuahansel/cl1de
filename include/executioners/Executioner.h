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
    const RunParameters & run_params);

protected:
  double computeElemSize() const;
  std::vector<double> computeElemPositions() const;
  std::vector<double> computeNodePositions() const;
  double getEndTime() const;
  double getTimeStepSizeAndUpdateTransientFlag(const double & dt_nominal, const double & t);
  std::vector<std::vector<double>> initializeRungeKuttaSolutionCoefs() const;
  std::vector<double> initializeRungeKuttaTimeStepCoefs() const;

  const Problem & _problem_base;
  const RunParameters & _run_params_base;

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
