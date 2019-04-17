#ifndef ExecutionerEuler1Phase_H
#define ExecutionerEuler1Phase_H

#include "Executioner.h"

class EOS1Phase;
class Function;
class ProblemEuler1Phase;
class RunParametersEuler1Phase;

class ExecutionerEuler1Phase : public Executioner
{
public:
  ExecutionerEuler1Phase(
    const ProblemEuler1Phase & problem,
    const RunParametersEuler1Phase & run_params);

  void run();

protected:
  std::vector<double> computeAreaNode() const;
  std::vector<double> computeAreaElem() const;

  const ProblemEuler1Phase & _problem;
  const RunParametersEuler1Phase & _run_params;

  const EOS1Phase & _eos;
  const Function & _A_fn;
  const Function & _r_ic_fn;
  const Function & _u_ic_fn;
  const Function & _p_ic_fn;

  const std::vector<double> _A_node;
  const std::vector<double> _A_elem;
};

#endif /* ExecutionerEuler1Phase_H */
