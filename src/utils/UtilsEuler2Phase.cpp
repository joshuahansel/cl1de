#include "UtilsEuler2Phase.h"

namespace Euler2Phase
{
  std::vector<double> computeConservativePhaseFlux(
    const std::vector<double> & W, const double & u_int, const double & p_int, const double & A)
  {
    const double & a = W[ia];
    const double & r = W[ir];
    const double & u = W[iu];
    const double & p = W[ip];
    const double & E = W[iE];

    std::vector<double> f(n_local_eq);

    f[0] = a * u_int * A;
    f[1] = a * r * u * A;
    f[2] = a * (r * u * u + p) * A - a * p_int * A;
    f[3] = a * u * (r * E + p) * A - a * p_int * u_int * A;

    return f;
  }

  std::vector<double> computeNonConservativePhaseFlux(
    const double & u_int, const double & p_int, const double & A)
  {
    std::vector<double> f(n_local_eq);

    f[0] = -u_int * A;
    f[2] = p_int * A;
    f[3] = p_int * u_int * A;

    return f;
  }

  std::vector<double> combineFluxes(
    const std::vector<double> & f_cons, const std::vector<double> & f_noncons, const double & a)
  {
    std::vector<double> f(f_cons.size());
    for (unsigned int i = 0; i < f_cons.size(); i++)
      f[i] = f_cons[i] + a * f_noncons[i];

    return f;
  }
}
