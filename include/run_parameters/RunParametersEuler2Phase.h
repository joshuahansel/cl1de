#ifndef RunParametersEuler2Phase_H
#define RunParametersEuler2Phase_H

#include "RunParameters.h"
#include "DoFHandlerEuler2Phase.h"

class EOS1Phase;
class FluxEuler2Phase;
class ReconstructorEuler2Phase;

class RunParametersEuler2Phase : public RunParameters
{
public:
  RunParametersEuler2Phase(int argc, char* argv[], const EOS1Phase & eos_liq, const EOS1Phase & eos_vap);

  const DoFHandlerEuler2Phase & getDoFHandler() const {return _dof_handler;}
  const FluxEuler2Phase & getFlux() const {return *_flux;}
  const ReconstructorEuler2Phase & getReconstructor() const {return *_reconstructor;}
  const bool & getUseFiniteVelocityRelaxation() const {return _use_finite_velocity_relaxation;}
  const bool & getUseFinitePressureRelaxation() const {return _use_finite_pressure_relaxation;}
  const bool & getUseInfiniteVelocityRelaxation() const {return _use_infinite_velocity_relaxation;}
  const bool & getUseInfinitePressureRelaxation() const {return _use_infinite_pressure_relaxation;}
  const unsigned int & getNumberOfSourceSubsteps() const {return _n_source_substeps;}

protected:
  const EOS1Phase & _eos_liq;
  const EOS1Phase & _eos_vap;
  const DoFHandlerEuler2Phase _dof_handler;

  std::shared_ptr<FluxEuler2Phase> _flux;
  std::shared_ptr<ReconstructorEuler2Phase> _reconstructor;

  bool _use_finite_velocity_relaxation;
  bool _use_finite_pressure_relaxation;

  bool _use_infinite_velocity_relaxation;
  bool _use_infinite_pressure_relaxation;

  unsigned int _n_source_substeps;
};

#endif /* RunParametersEuler2Phase_H */
