#ifndef BCEuler2Phase_H
#define BCEuler2Phase_H

#include "BC.h"

#include <vector>

class BCEuler2Phase : public BC
{
public:
  BCEuler2Phase(bool is_left);

  virtual std::vector<double> computeFlux(
    const std::vector<double> & U, const double & A) const = 0;
};

#endif /* BCEuler2Phase_H */
