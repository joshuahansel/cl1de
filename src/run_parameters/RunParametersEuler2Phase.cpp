#include "RunParametersEuler2Phase.h"
#include "EOS1Phase.h"
#include "FluxEuler2PhaseHLLC.h"
#include "ReconstructorEuler2PhaseGodunov.h"
#include "ReconstructorEuler2PhaseSlopePUTFull.h"
#include "ReconstructorEuler2PhaseSlopePUTMC.h"
#include "ReconstructorEuler2PhaseSlopePUTMinmod.h"
#include "ReconstructorEuler2PhaseSlopePUTSuperbee.h"
#include "utils.h"

RunParametersEuler2Phase::RunParametersEuler2Phase(
  int argc, char* argv[], const EOS1Phase & eos_liq, const EOS1Phase & eos_vap)
  : RunParameters(argc, argv),
    _eos_liq(eos_liq),
    _eos_vap(eos_vap),
    _dof_handler(DoFHandlerEuler2Phase(_n_elem))
{
  // flux
  const std::string flux_option = getStringParameter("flux");
  if (flux_option == "hllc")
    _flux = std::make_shared<FluxEuler2PhaseHLLC>(_eos_liq, _eos_vap);
  else
    throwInvalidStringParameterValueError("flux", flux_option);

  // reconstructor
  const std::string slope_reconstruction = getStringParameter("slope_reconstruction");
  if (slope_reconstruction == "none")
    _reconstructor = std::make_shared<ReconstructorEuler2PhaseGodunov>(_dof_handler);
  else if (slope_reconstruction == "full")
    _reconstructor = std::make_shared<ReconstructorEuler2PhaseSlopePUTFull>(_dof_handler, _eos_liq, _eos_vap);
  else if (slope_reconstruction == "mc")
    _reconstructor = std::make_shared<ReconstructorEuler2PhaseSlopePUTMC>(_dof_handler, _eos_liq, _eos_vap);
  else if (slope_reconstruction == "minmod")
    _reconstructor = std::make_shared<ReconstructorEuler2PhaseSlopePUTMinmod>(_dof_handler, _eos_liq, _eos_vap);
  else if (slope_reconstruction == "superbee")
    _reconstructor = std::make_shared<ReconstructorEuler2PhaseSlopePUTSuperbee>(_dof_handler, _eos_liq, _eos_vap);
  else
    throwInvalidStringParameterValueError("slope_reconstruction", slope_reconstruction);

  // finite relaxation
  _use_finite_velocity_relaxation = getBoolParameter("use_finite_velocity_relaxation");
  _use_finite_pressure_relaxation = getBoolParameter("use_finite_pressure_relaxation");

  // infinite relaxation
  _use_infinite_velocity_relaxation = getBoolParameter("use_infinite_velocity_relaxation");
  _use_infinite_pressure_relaxation = getBoolParameter("use_infinite_pressure_relaxation");

  // number of source substeps
  _n_source_substeps = getIntParameter("n_source_substeps");

  checkAllParametersUsed();
}
