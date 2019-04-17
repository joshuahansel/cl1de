#ifndef BCEuler1Phase_H
#define BCEuler1Phase_H

#include "BC.h"

#include <vector>

class BCEuler1Phase : public BC
{
public:
  BCEuler1Phase(bool is_left);

  virtual std::vector<double> computeFlux(
    const std::vector<double> & U, const double & A) const = 0;
};

#endif /* BCEuler1Phase_H */
