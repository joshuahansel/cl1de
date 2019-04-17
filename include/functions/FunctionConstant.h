#ifndef FunctionConstant_H
#define FunctionConstant_H

#include "Function.h"

class FunctionConstant : public Function
{
public:
  FunctionConstant(const double & value);

  virtual double value(const double & x, const double & t) const override;

protected:
  const double _value;
};

#endif /* FunctionConstant_H */
