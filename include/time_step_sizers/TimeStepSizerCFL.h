#ifndef TimeStepSizerCFL_H
#define TimeStepSizerCFL_H

#include "TimeStepSizer.h"

class TimeStepSizerCFL : public TimeStepSizer
{
public:
  TimeStepSizerCFL(const double & cfl);

  virtual double computeTimeStepSize(const double & max_wave_speed, const double & dx_min) const override;

protected:
  const double _cfl;
};

#endif /* TimeStepSizerCFL_H */
