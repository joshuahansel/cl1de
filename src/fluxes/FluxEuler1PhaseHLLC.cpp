#include "FluxEuler1PhaseHLLC.h"
#include "UtilsEuler1Phase.h"
#include "UtilsNumerics.h"
#include "EOS1Phase.h"

#include <cmath>

FluxEuler1PhaseHLLC::FluxEuler1PhaseHLLC(const EOS1Phase & eos)
  : FluxEuler1Phase(eos)
{
}

std::vector<double> FluxEuler1PhaseHLLC::computeFlux(
  const std::vector<double> & U_L,
  const std::vector<double> & U_R,
  const double & A) const
{
  const double & rA_L = U_L[0];
  const double & ruA_L = U_L[1];
  const double & rEA_L = U_L[2];
  const double r_L = rA_L / A;
  const double u_L = ruA_L / rA_L;
  const double E_L = rEA_L / rA_L;
  const double e_L = E_L - 0.5 * u_L * u_L;
  const double p_L = _eos.p_from_r_e(r_L, e_L);
  const double c_L = _eos.c_from_r_e(r_L, e_L);
  const double H_L = E_L + p_L / r_L;

  const double & rA_R = U_R[0];
  const double & ruA_R = U_R[1];
  const double & rEA_R = U_R[2];
  const double r_R = rA_R / A;
  const double u_R = ruA_R / rA_R;
  const double E_R = rEA_R / rA_R;
  const double e_R = E_R - 0.5 * u_R * u_R;
  const double p_R = _eos.p_from_r_e(r_R, e_R);
  const double c_R = _eos.c_from_r_e(r_R, e_R);
  const double H_R = E_R + p_R / r_R;

  // compute wave speeds
  // const double r_roe = std::sqrt(r_L * r_R);
  // const double sqrt_r_L = std::sqrt(r_L);
  // const double sqrt_r_R = std::sqrt(r_R);
  // const double u_roe = Numerics::weightedAverage(u_L, u_R, sqrt_r_L, sqrt_r_R);
  // const double H_roe = Numerics::weightedAverage(H_L, H_R, sqrt_r_L, sqrt_r_R);
  // const double h_roe = H_roe - 0.5 * u_roe * u_roe;
  // const double e_roe = _eos.e_from_r_h(r_roe, h_roe);
  // const double c_roe = _eos.c_from_r_e(r_roe, e_roe);
  // const double lambda1_L = u_L - c_L;
  // const double lambda3_R = u_R + c_R;
  // const double lambda1_roe = u_roe - c_roe;
  // const double lambda3_roe = u_roe + c_roe;
  // const double S_L = std::min(lambda1_L, lambda1_roe);
  // const double S_R = std::max(lambda3_R, lambda3_roe);
  const double S_L = std::min(u_L - c_L, u_R - c_R);
  const double S_R = std::max(u_L + c_L, u_R + c_R);

  const double S_M = (r_R * u_R * (S_R - u_R) - r_L * u_L * (S_L - u_L) + p_L - p_R)
    / (r_R * (S_R - u_R) - r_L * (S_L - u_L));
  const double p_star = r_L * (u_L - S_L) * (u_L - S_M) + p_L;

  if (S_L > 0)
  {
    return Euler1Phase::computeFlux(U_L, A, _eos);
  }
  else if (S_M > 0)
  {
    const double factor = A / (S_L - S_M);
    const double rA = (S_L - u_L) * r_L * factor;
    const double ruA = ((S_L - u_L) * r_L * u_L + (p_star - p_L)) * factor;
    const double rEA = ((S_L - u_L) * r_L * E_L - p_L * u_L + p_star * S_M) * factor;

    std::vector<double> f(3, 0.0);
    f[0] = rA * S_M;
    f[1] = ruA * S_M + p_star * A;
    f[2] = rEA * S_M + p_star * S_M * A;
    return f;
  }
  else if (S_R > 0)
  {
    const double factor = A / (S_R - S_M);
    const double rA = (S_R - u_R) * r_R * factor;
    const double ruA = ((S_R - u_R) * r_R * u_R + (p_star - p_R)) * factor;
    const double rEA = ((S_R - u_R) * r_R * E_R - p_R * u_R + p_star * S_M) * factor;

    std::vector<double> f(3, 0.0);
    f[0] = rA * S_M;
    f[1] = ruA * S_M + p_star * A;
    f[2] = rEA * S_M + p_star * S_M * A;
    return f;
  }
  else
  {
    return Euler1Phase::computeFlux(U_R, A, _eos);
  }
}
