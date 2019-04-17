#ifndef TimeStepSizer_H
#define TimeStepSizer_H

class TimeStepSizer
{
public:
  TimeStepSizer();

  virtual double computeTimeStepSize(const double & max_wave_speed, const double & dx_min) const = 0;
};

#endif /* TimeStepSizer_H */
