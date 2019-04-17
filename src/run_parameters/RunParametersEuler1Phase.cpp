#include "RunParametersEuler1Phase.h"
#include "EOS1Phase.h"
#include "FluxEuler1PhaseCentered.h"
#include "FluxEuler1PhaseHLLC.h"
#include "ReconstructorEuler1PhaseGodunov.h"
#include "ReconstructorEuler1PhaseSlopeConservativeMC.h"
#include "ReconstructorEuler1PhaseSlopeConservativeMinmod.h"
#include "ReconstructorEuler1PhaseSlopeConservativeSuperbee.h"
#include "ReconstructorEuler1PhaseSlopePUTMC.h"
#include "ReconstructorEuler1PhaseSlopePUTMinmod.h"
#include "ReconstructorEuler1PhaseSlopePUTSuperbee.h"
#include "ReconstructorEuler1PhaseSlopeRUPMC.h"
#include "ReconstructorEuler1PhaseSlopeRUPMinmod.h"
#include "ReconstructorEuler1PhaseSlopeRUPSuperbee.h"
#include "utils.h"

RunParametersEuler1Phase::RunParametersEuler1Phase(int argc, char* argv[], const EOS1Phase & eos)
  : RunParameters(argc, argv),
    _eos(eos)
{
  // flux
  const std::string flux_option = getStringParameter("flux");
  if (flux_option == "centered")
    _flux = std::make_shared<FluxEuler1PhaseCentered>(_eos);
  else if (flux_option == "hllc")
    _flux = std::make_shared<FluxEuler1PhaseHLLC>(_eos);
  else
    throwInvalidStringParameterValueError("flux", flux_option);

  // reconstructor
  const std::string slope_reconstruction = getStringParameter("slope_reconstruction");
  if (slope_reconstruction == "none")
    _reconstructor = std::make_shared<ReconstructorEuler1PhaseGodunov>();
  else if (slope_reconstruction == "minmod")
  {
    const std::string reconstruction_var_set = getStringParameter("reconstruction_var_set");
    if (reconstruction_var_set == "conservative")
      _reconstructor = std::make_shared<ReconstructorEuler1PhaseSlopeConservativeMinmod>(_eos);
    else if (reconstruction_var_set == "puT")
      _reconstructor = std::make_shared<ReconstructorEuler1PhaseSlopePUTMinmod>(_eos);
    else if (reconstruction_var_set == "rup")
      _reconstructor = std::make_shared<ReconstructorEuler1PhaseSlopeRUPMinmod>(_eos);
    else
      throwInvalidStringParameterValueError("reconstruction_var_set", reconstruction_var_set);
  }
  else if (slope_reconstruction == "mc")
  {
    const std::string reconstruction_var_set = getStringParameter("reconstruction_var_set");
    if (reconstruction_var_set == "conservative")
      _reconstructor = std::make_shared<ReconstructorEuler1PhaseSlopeConservativeMC>(_eos);
    else if (reconstruction_var_set == "puT")
      _reconstructor = std::make_shared<ReconstructorEuler1PhaseSlopePUTMC>(_eos);
    else if (reconstruction_var_set == "rup")
      _reconstructor = std::make_shared<ReconstructorEuler1PhaseSlopeRUPMC>(_eos);
    else
      throwInvalidStringParameterValueError("reconstruction_var_set", reconstruction_var_set);
  }
  else if (slope_reconstruction == "superbee")
  {
    const std::string reconstruction_var_set = getStringParameter("reconstruction_var_set");
    if (reconstruction_var_set == "conservative")
      _reconstructor = std::make_shared<ReconstructorEuler1PhaseSlopeConservativeSuperbee>(_eos);
    else if (reconstruction_var_set == "puT")
      _reconstructor = std::make_shared<ReconstructorEuler1PhaseSlopePUTSuperbee>(_eos);
    else if (reconstruction_var_set == "rup")
      _reconstructor = std::make_shared<ReconstructorEuler1PhaseSlopeRUPSuperbee>(_eos);
    else
      throwInvalidStringParameterValueError("reconstruction_var_set", reconstruction_var_set);
  }
  else
    throwInvalidStringParameterValueError("slope_reconstruction", slope_reconstruction);

  checkAllParametersUsed();
}
