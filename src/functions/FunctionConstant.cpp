#include "FunctionConstant.h"

FunctionConstant::FunctionConstant(const double & value) : Function(), _value(value)
{
}

double FunctionConstant::value(const double & /*x*/, const double & /*t*/) const
{
  return _value;
}
