#include "TimeStepSizerConstant.h"

TimeStepSizerConstant::TimeStepSizerConstant(const double & dt)
  : TimeStepSizer(), _dt(dt)
{
}

double TimeStepSizerConstant::computeTimeStepSize(
  const double & /*max_wave_speed*/, const double & /*dx_min*/) const
{
  return _dt;
}
