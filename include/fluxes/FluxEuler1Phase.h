#ifndef FluxEuler1Phase_H
#define FluxEuler1Phase_H

#include <vector>

class EOS1Phase;

class FluxEuler1Phase
{
public:
  FluxEuler1Phase(const EOS1Phase & eos);

  virtual std::vector<double> computeFlux(
    const std::vector<double> & U_L,
    const std::vector<double> & U_R,
    const double & A) const = 0;

protected:
  const EOS1Phase & _eos;
};

#endif /* FluxEuler1Phase_H */
