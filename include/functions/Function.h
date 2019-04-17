#ifndef Function_H
#define Function_H

class Function
{
public:
  Function();

  virtual double value(const double & x, const double & t) const = 0;
};

#endif /* Function_H */
