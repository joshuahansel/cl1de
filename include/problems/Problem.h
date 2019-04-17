#ifndef Problem_H
#define Problem_H

class Problem
{
public:
  Problem();

  void setDomain(const double & x_min, const double & x_max);
  void setDefaultEndTime(const double & t_end);

  const double & getXMin() const {return _x_min;}
  const double & getXMax() const {return _x_max;}
  const double & getDefaultEndTime() const {return _t_end;}

  bool hasDefaultEndTime() const {return _t_end_set;}

protected:
  double _x_min;
  double _x_max;
  double _t_end;

  bool _t_end_set;
};

#endif /* Problem_H */
