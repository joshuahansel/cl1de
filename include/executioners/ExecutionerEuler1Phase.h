#ifndef ExecutionerEuler1Phase_H
#define ExecutionerEuler1Phase_H

#include "Executioner.h"

class BCEuler1Phase;
class EOS1Phase;
class FluxEuler1Phase;
class Function;
class ICsEuler1Phase;
class ProblemEuler1Phase;
class ReconstructorEuler1Phase;
class RunParametersEuler1Phase;
class DoFHandlerEuler1Phase;

class ExecutionerEuler1Phase : public Executioner
{
public:
  ExecutionerEuler1Phase(
    const ProblemEuler1Phase & problem,
    const RunParametersEuler1Phase & run_params);

protected:
  virtual void initializeSolution(std::vector<double> & U) const override;
  virtual double computeMaxWaveSpeed(const std::vector<double> & U) const override;
  virtual void computeSteadyStateResidual(
    const std::vector<double> & U, std::vector<double> & ss_rhs) const override;
  virtual void outputSolution(const std::vector<double> & U) const override;

  std::vector<double> computeAreaNode() const;
  std::vector<double> computeAreaElem() const;

  const ProblemEuler1Phase & _problem;
  const RunParametersEuler1Phase & _run_params;

  const DoFHandlerEuler1Phase & _dof_handler;
  const EOS1Phase & _eos;
  const Function & _A_fn;
  const ICsEuler1Phase & _ics;
  const BCEuler1Phase & _bc_left;
  const BCEuler1Phase & _bc_right;
  const FluxEuler1Phase & _flux;
  const ReconstructorEuler1Phase & _reconstructor;

  const std::vector<double> _A_node;
  const std::vector<double> _A_elem;
};

#endif /* ExecutionerEuler1Phase_H */
