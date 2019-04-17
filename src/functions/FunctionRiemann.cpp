#include "FunctionRiemann.h"

FunctionRiemann::FunctionRiemann(
  const double & left_value,
  const double & right_value,
  const double & x_interface)
  : Function(),
    _left_value(left_value),
    _right_value(right_value),
    _x_interface(x_interface)
{
}

double FunctionRiemann::value(const double & x, const double & /*t*/) const
{
  if (x < _x_interface)
    return _left_value;
  else
    return _right_value;
}
