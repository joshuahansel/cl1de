#ifndef UtilsEuler2Phase_H
#define UtilsEuler2Phase_H

#include <vector>

class EOS1Phase;

namespace Euler2Phase
{
  const unsigned int n_eq = 7;
  const unsigned int n_local_eq = 4;

  const unsigned int n_local_prim_var = 5;
  const unsigned int ia = 0;
  const unsigned int ir = 1;
  const unsigned int iu = 2;
  const unsigned int ip = 3;
  const unsigned int iE = 4;

  std::vector<double> computeConservativePhaseFlux(
    const std::vector<double> & W, const double & u_int, const double & p_int, const double & A);

  std::vector<double> computeNonConservativePhaseFlux(
    const double & u_int, const double & p_int, const double & A);

  std::vector<double> combineFluxes(
    const std::vector<double> & f_cons, const std::vector<double> & f_noncons, const double & a);
}

#endif /* UtilsEuler2Phase_H */
