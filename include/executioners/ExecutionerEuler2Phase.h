#ifndef ExecutionerEuler2Phase_H
#define ExecutionerEuler2Phase_H

#include "Executioner.h"

class BCEuler2Phase;
class EOS1Phase;
class FluxEuler2Phase;
class Function;
class ICsEuler2Phase;
class ProblemEuler2Phase;
class ReconstructorEuler2Phase;
class RunParametersEuler2Phase;
class DoFHandlerEuler2Phase;

class ExecutionerEuler2Phase : public Executioner
{
public:
  ExecutionerEuler2Phase(
    const ProblemEuler2Phase & problem,
    const RunParametersEuler2Phase & run_params);

protected:
  virtual void initializeSolution(std::vector<double> & U) const override;
  virtual double computeMaxWaveSpeed(const std::vector<double> & U) const override;
  virtual void computeSteadyStateResidual(
    const std::vector<double> & U, std::vector<double> & ss_rhs) const override;
  virtual void performPostStep(std::vector<double> & U, const double & dt) const override;
  virtual void outputSolution(const std::vector<double> & U) const override;

  std::vector<double> computeAreaNode() const;
  std::vector<double> computeAreaElem() const;

  const ProblemEuler2Phase & _problem;
  const RunParametersEuler2Phase & _run_params;

  const DoFHandlerEuler2Phase & _dof_handler;
  const EOS1Phase & _eos_liq;
  const EOS1Phase & _eos_vap;
  const Function & _A_fn;
  const ICsEuler2Phase & _ics;
  const BCEuler2Phase & _bc_left;
  const BCEuler2Phase & _bc_right;
  const double & _gravity;
  const double & _a_int;
  const FluxEuler2Phase & _flux;
  const ReconstructorEuler2Phase & _reconstructor;

  const bool & _use_finite_velocity_relaxation;
  const bool & _use_finite_pressure_relaxation;

  const bool & _use_infinite_velocity_relaxation;
  const bool & _use_infinite_pressure_relaxation;

  const unsigned int & _n_source_substeps;

  const std::vector<double> _A_node;
  const std::vector<double> _A_elem;
};

#endif /* ExecutionerEuler2Phase_H */
