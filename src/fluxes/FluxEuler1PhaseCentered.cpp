#include "FluxEuler1PhaseCentered.h"
#include "UtilsEuler1Phase.h"

FluxEuler1PhaseCentered::FluxEuler1PhaseCentered(const EOS1Phase & eos)
  : FluxEuler1Phase(eos)
{
}

std::vector<double> FluxEuler1PhaseCentered::computeFlux(
  const std::vector<double> & U_L,
  const std::vector<double> & U_R,
  const double & A) const
{
  const std::vector<double> f_L = Euler1Phase::computeFlux(U_L, A, _eos);
  const std::vector<double> f_R = Euler1Phase::computeFlux(U_R, A, _eos);

  std::vector<double> f(3, 0.0);
  for (unsigned int i = 0; i < 3; i++)
    f[i] = 0.5 * (f_L[i] + f_R[i]);

  return f;
}
