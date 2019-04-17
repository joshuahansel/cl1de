#include "UtilsEuler1Phase.h"
#include "EOS1Phase.h"

namespace Euler1Phase
{
  std::vector<double> computeFlux(
    const std::vector<double> & U, const double & A, const EOS1Phase & eos)
  {
    const double r = U[0] / A;
    const double u = U[1] / U[0];
    const double E = U[2] / U[0];
    const double e = E - 0.5 * u * u;

    const double p = eos.p_from_r_e(r, e);

    std::vector<double> f(3, 0.0);
    f[0] = U[1];
    f[1] = (r * u * u + p) * A;
    f[2] = u * (r * E + p) * A;

    return f;
  }
}
