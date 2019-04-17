#ifndef TimeStepSizerConstant_H
#define TimeStepSizerConstant_H

#include "TimeStepSizer.h"

class TimeStepSizerConstant : public TimeStepSizer
{
public:
  TimeStepSizerConstant(const double & dt);

  virtual double computeTimeStepSize(const double & max_wave_speed, const double & dx_min) const override;

protected:
  const double _dt;
};

#endif /* TimeStepSizerConstant_H */
