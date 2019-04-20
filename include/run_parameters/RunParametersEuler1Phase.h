#ifndef RunParametersEuler1Phase_H
#define RunParametersEuler1Phase_H

#include "RunParameters.h"
#include "DoFHandlerEuler1Phase.h"

class EOS1Phase;
class FluxEuler1Phase;
class ReconstructorEuler1Phase;

class RunParametersEuler1Phase : public RunParameters
{
public:
  RunParametersEuler1Phase(int argc, char* argv[], const EOS1Phase & eos);

  const DoFHandlerEuler1Phase & getDoFHandler() const {return _dof_handler;}
  const FluxEuler1Phase & getFlux() const {return *_flux;}
  const ReconstructorEuler1Phase & getReconstructor() const {return *_reconstructor;}

protected:
  const EOS1Phase & _eos;
  const DoFHandlerEuler1Phase _dof_handler;

  std::shared_ptr<FluxEuler1Phase> _flux;
  std::shared_ptr<ReconstructorEuler1Phase> _reconstructor;
};

#endif /* RunParametersEuler1Phase_H */
