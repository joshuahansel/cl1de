#ifndef FunctionRiemann_H
#define FunctionRiemann_H

#include "Function.h"

class FunctionRiemann : public Function
{
public:
  FunctionRiemann(const double & left_value, const double & right_value, const double & x_interface);

  virtual double value(const double & x, const double & t) const override;

protected:
  const double _left_value;
  const double _right_value;
  const double _x_interface;
};

#endif /* FunctionRiemann_H */
