#include "TimeStepSizerCFL.h"

TimeStepSizerCFL::TimeStepSizerCFL(const double & cfl)
  : TimeStepSizer(), _cfl(cfl)
{
}

double TimeStepSizerCFL::computeTimeStepSize(
  const double & max_wave_speed, const double & dx_min) const
{
  return _cfl * dx_min / max_wave_speed;
}
