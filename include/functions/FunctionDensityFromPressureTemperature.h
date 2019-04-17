#ifndef FunctionDensityFromPressureTemperature_H
#define FunctionDensityFromPressureTemperature_H

#include "Function.h"

class EOS1Phase;

class FunctionDensityFromPressureTemperature : public Function
{
public:
  FunctionDensityFromPressureTemperature(
    const EOS1Phase & eos,
    const Function & pressure_fn,
    const Function & temperature_fn);

  virtual double value(const double & x, const double & t) const override;

protected:
  const EOS1Phase & _eos;
  const Function & _pressure_fn;
  const Function & _temperature_fn;
};

#endif /* FunctionDensityFromPressureTemperature_H */
