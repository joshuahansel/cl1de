#ifndef FunctionCosineHump_H
#define FunctionCosineHump_H

#include "Function.h"

class FunctionCosineHump : public Function
{
public:
  FunctionCosineHump(
    const double & value_begin,
    const double & value_center,
    const double & hump_center,
    const double & hump_width);

  virtual double value(const double & x, const double & t) const override;

protected:
  const double _value_begin;
  const double _value_center;
  const double _hump_center;
  const double _hump_width;

  const double _hump_begin;
  const double _hump_end;
  const double _cosine_amplitude;
  const double _value_mid;
};

#endif /* FunctionCosineHump_H */
