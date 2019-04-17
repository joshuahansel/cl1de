#include "FunctionDensityFromPressureTemperature.h"
#include "EOS1Phase.h"

FunctionDensityFromPressureTemperature::FunctionDensityFromPressureTemperature(
  const EOS1Phase & eos,
  const Function & pressure_fn,
  const Function & temperature_fn)
  : Function(),
    _eos(eos),
    _pressure_fn(pressure_fn),
    _temperature_fn(temperature_fn)
{
}

double FunctionDensityFromPressureTemperature::value(const double & x, const double & t) const
{
  return _eos.r_from_p_T(_pressure_fn.value(x, t), _temperature_fn.value(x, t));
}
