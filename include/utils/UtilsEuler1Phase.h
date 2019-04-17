#ifndef UtilsEuler1Phase_H
#define UtilsEuler1Phase_H

#include <vector>

class EOS1Phase;

namespace Euler1Phase
{
  std::vector<double> computeFlux(
    const std::vector<double> & U, const double & A, const EOS1Phase & eos);
}

#endif /* UtilsEuler1Phase_H */
