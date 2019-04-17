#include "BCEuler1PhaseFree.h"
#include "UtilsEuler1Phase.h"

BCEuler1PhaseFree::BCEuler1PhaseFree(bool is_left, const EOS1Phase & eos)
  : BCEuler1Phase(is_left),
    _eos(eos)
{
}

std::vector<double> BCEuler1PhaseFree::computeFlux(
  const std::vector<double> & U, const double & A) const
{
  return Euler1Phase::computeFlux(U, A, _eos);
}
