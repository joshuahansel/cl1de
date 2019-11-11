#include "FunctionCosineHump.h"
#include "PhysicsConstants.h"

#include <cmath>

FunctionCosineHump::FunctionCosineHump(
  const double & value_begin,
  const double & value_center,
  const double & hump_center,
  const double & hump_width)
  : Function(),

    _value_begin(value_begin),
    _value_center(value_center),
    _hump_center(hump_center),
    _hump_width(hump_width),

    _hump_begin(_hump_center - 0.5 * _hump_width),
    _hump_end(_hump_center + 0.5 * _hump_width),
    _cosine_amplitude(0.5 * (_value_center - _value_begin)),
    _value_mid(0.5 * (_value_center + _value_begin))
{
}

double FunctionCosineHump::value(const double & x, const double & /*t*/) const
{
  if (x <= _hump_begin)
    return _value_begin;
  else if (x >= _hump_end)
    return _value_begin;
  else
    return _value_mid - _cosine_amplitude * std::cos(
      2 * PhysicsConstants::pi / _hump_width * (x - _hump_center + 0.5 * _hump_width));
}
