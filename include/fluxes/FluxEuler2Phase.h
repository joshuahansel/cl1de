#ifndef FluxEuler2Phase_H
#define FluxEuler2Phase_H

#include <vector>

class EOS1Phase;

class FluxEuler2Phase
{
public:
  FluxEuler2Phase(const EOS1Phase & eos_liq, const EOS1Phase & eos_vap);

  virtual void computeFlux(
    const std::vector<double> & U_L,
    const std::vector<double> & U_R,
    const double & A,
    std::vector<double> & f_L,
    std::vector<double> & f_R) const = 0;

protected:
  const EOS1Phase & _eos_liq;
  const EOS1Phase & _eos_vap;
  const std::vector<const EOS1Phase *> _eos;
};

#endif /* FluxEuler2Phase_H */
