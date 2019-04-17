#include "Problem.h"

Problem::Problem() : _t_end_set(false)
{
}

void Problem::setDomain(const double & x_min, const double & x_max)
{
  _x_min = x_min;
  _x_max = x_max;
}

void Problem::setDefaultEndTime(const double & t_end)
{
  _t_end = t_end;

  _t_end_set = true;
}
