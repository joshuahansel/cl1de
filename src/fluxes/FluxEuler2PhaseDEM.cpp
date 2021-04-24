#include "FluxEuler2PhaseDEM.h"
#include "UtilsEuler1Phase.h"
#include "EOS1Phase.h"

FluxEuler2PhaseDEM::FluxEuler2PhaseDEM(const EOS1Phase & eos_liq, const EOS1Phase & eos_vap)
  : FluxEuler2Phase(eos_liq, eos_vap)
{
}

void FluxEuler2PhaseDEM::computeFlux(
  const std::vector<double> & U_L,
  const std::vector<double> & U_R,
  const double & A,
  std::vector<double> & f_L,
  std::vector<double> & f_R) const
{
  const std::vector<std::vector<double>> U {U_L, U_R};

  std::vector<std::vector<double>> a(_n_phases, std::vector<double>(_n_sides));
  std::vector<std::vector<double>> U_liq(_n_sides, std::vector<double>(3));
  std::vector<std::vector<double>> U_vap(_n_sides, std::vector<double>(3));
  for (unsigned int i = 0; i < _n_sides; i++)
  {
    a[liq][i] = U[i][0] / A;
    a[vap][i] = 1 - a[liq][i];

    U_liq[i][0] = U[i][1] / a[liq][i];
    U_liq[i][1] = U[i][2] / a[liq][i];
    U_liq[i][2] = U[i][3] / a[liq][i];

    U_vap[i][0] = U[i][4] / a[vap][i];
    U_vap[i][1] = U[i][5] / a[vap][i];
    U_vap[i][2] = U[i][6] / a[vap][i];
  }

  const double S_liq_liq = std::min(a[liq][L], a[liq][R]);
  const double S_liq_vap = std::max(a[liq][L] - a[liq][R], 0.0);
  const double S_vap_liq = std::max(a[vap][L] - a[vap][R], 0.0);
  const double S_vap_vap = std::min(a[vap][L], a[vap][R]);

  if (S_liq_liq > 0)
  {
    double u_int, p_int;
    const auto f_conv = computeConvectiveFlux1Phase(U_liq[L], U_liq[R], _eos_liq, _eos_liq, A, u_int, p_int);
    apply1PhaseFlux(f_conv, S_liq_liq, liq, f_L);
    apply1PhaseFlux(f_conv, S_liq_liq, liq, f_R);
  }
  if (S_liq_vap > 0)
  {
    double u_int, p_int;
    const auto f_conv = computeConvectiveFlux1Phase(U_liq[L], U_vap[R], _eos_liq, _eos_vap, A, u_int, p_int);

    double x_liq, x_jump_liq_L, x_jump_liq_R, x_jump_vap_L, x_jump_vap_R;
    if (u_int > 0)
    {
      x_liq = 1;
      x_jump_liq_L = 0;
      x_jump_liq_R = -1;
      x_jump_vap_L = 0;
      x_jump_vap_R = 1;
    }
    else
    {
      x_liq = 0;
      x_jump_liq_L = -1;
      x_jump_liq_R = 0;
      x_jump_vap_L = 1;
      x_jump_vap_R = 0;
    }
    const double x_vap = 1 - x_liq;

    apply1PhaseFlux(f_conv, x_liq * S_liq_vap, liq, f_L);
    apply1PhaseFlux(f_conv, x_liq * S_liq_vap, liq, f_R);
    apply1PhaseFlux(f_conv, x_vap * S_liq_vap, vap, f_L);
    apply1PhaseFlux(f_conv, x_vap * S_liq_vap, vap, f_R);

    std::vector<double> f_lag_L(7);
    f_lag_L[0] = -u_int * x_jump_liq_L;
    f_lag_L[2] = p_int * x_jump_liq_L;
    f_lag_L[3] = u_int * p_int * x_jump_liq_L;
    f_lag_L[5] = p_int * x_jump_vap_L;
    f_lag_L[6] = u_int * p_int * x_jump_vap_L;
    apply2PhaseFlux(f_lag_L, -S_liq_vap * A, f_L);

    std::vector<double> f_lag_R(7);
    f_lag_R[0] = -u_int * x_jump_liq_R;
    f_lag_R[2] = p_int * x_jump_liq_R;
    f_lag_R[3] = u_int * p_int * x_jump_liq_R;
    f_lag_R[5] = p_int * x_jump_vap_R;
    f_lag_R[6] = u_int * p_int * x_jump_vap_R;
    apply2PhaseFlux(f_lag_R, S_liq_vap * A, f_R);
  }
  if (S_vap_liq > 0)
  {
    double u_int, p_int;
    const auto f_conv = computeConvectiveFlux1Phase(U_vap[L], U_liq[R], _eos_vap, _eos_liq, A, u_int, p_int);

    double x_liq, x_jump_liq_L, x_jump_liq_R, x_jump_vap_L, x_jump_vap_R;
    if (u_int > 0)
    {
      x_liq = 0;
      x_jump_liq_L = 0;
      x_jump_liq_R = 1;
      x_jump_vap_L = 0;
      x_jump_vap_R = -1;
    }
    else
    {
      x_liq = 1;
      x_jump_liq_L = 1;
      x_jump_liq_R = 0;
      x_jump_vap_L = -1;
      x_jump_vap_R = 0;
    }
    const double x_vap = 1 - x_liq;

    apply1PhaseFlux(f_conv, x_liq * S_vap_liq, liq, f_L);
    apply1PhaseFlux(f_conv, x_liq * S_vap_liq, liq, f_R);
    apply1PhaseFlux(f_conv, x_vap * S_vap_liq, vap, f_L);
    apply1PhaseFlux(f_conv, x_vap * S_vap_liq, vap, f_R);

    std::vector<double> f_lag_L(7);
    f_lag_L[0] = -u_int * x_jump_liq_L;
    f_lag_L[2] = p_int * x_jump_liq_L;
    f_lag_L[3] = u_int * p_int * x_jump_liq_L;
    f_lag_L[5] = p_int * x_jump_vap_L;
    f_lag_L[6] = u_int * p_int * x_jump_vap_L;
    apply2PhaseFlux(f_lag_L, -S_vap_liq * A, f_L);

    std::vector<double> f_lag_R(7);
    f_lag_R[0] = -u_int * x_jump_liq_R;
    f_lag_R[2] = p_int * x_jump_liq_R;
    f_lag_R[3] = u_int * p_int * x_jump_liq_R;
    f_lag_R[5] = p_int * x_jump_vap_R;
    f_lag_R[6] = u_int * p_int * x_jump_vap_R;
    apply2PhaseFlux(f_lag_R, S_vap_liq * A, f_R);
  }
  if (S_vap_vap > 0)
  {
    double u_int, p_int;
    const auto f_conv = computeConvectiveFlux1Phase(U_vap[L], U_vap[R], _eos_vap, _eos_vap, A, u_int, p_int);
    apply1PhaseFlux(f_conv, S_vap_vap, vap, f_L);
    apply1PhaseFlux(f_conv, S_vap_vap, vap, f_R);
  }
}

std::vector<double> FluxEuler2PhaseDEM::computeConvectiveFlux1Phase(
  const std::vector<double> & U_L,
  const std::vector<double> & U_R,
  const EOS1Phase & eos_L,
  const EOS1Phase & eos_R,
  const double & A,
  double & u_int,
  double & p_int) const
{
  const double & rA_L = U_L[0];
  const double & ruA_L = U_L[1];
  const double & rEA_L = U_L[2];
  const double r_L = rA_L / A;
  const double u_L = ruA_L / rA_L;
  const double E_L = rEA_L / rA_L;
  const double e_L = E_L - 0.5 * u_L * u_L;
  const double p_L = eos_L.p_from_r_e(r_L, e_L);
  const double c_L = eos_L.c_from_r_e(r_L, e_L);
  const double H_L = E_L + p_L / r_L;

  const double & rA_R = U_R[0];
  const double & ruA_R = U_R[1];
  const double & rEA_R = U_R[2];
  const double r_R = rA_R / A;
  const double u_R = ruA_R / rA_R;
  const double E_R = rEA_R / rA_R;
  const double e_R = E_R - 0.5 * u_R * u_R;
  const double p_R = eos_R.p_from_r_e(r_R, e_R);
  const double c_R = eos_R.c_from_r_e(r_R, e_R);
  const double H_R = E_R + p_R / r_R;

  // compute wave speeds
  const double S_L = std::min(u_L - c_L, u_R - c_R);
  const double S_R = std::max(u_L + c_L, u_R + c_R);

  u_int = (r_R * u_R * (S_R - u_R) - r_L * u_L * (S_L - u_L) + p_L - p_R)
    / (r_R * (S_R - u_R) - r_L * (S_L - u_L));
  p_int = r_L * (u_L - S_L) * (u_L - u_int) + p_L;

  if (S_L > 0)
  {
    return Euler1Phase::computeFlux(U_L, A, eos_L);
  }
  else if (u_int > 0)
  {
    const double factor = A / (S_L - u_int);
    const double rA = (S_L - u_L) * r_L * factor;
    const double ruA = ((S_L - u_L) * r_L * u_L + (p_int - p_L)) * factor;
    const double rEA = ((S_L - u_L) * r_L * E_L - p_L * u_L + p_int * u_int) * factor;

    std::vector<double> f(3, 0.0);
    f[0] = rA * u_int;
    f[1] = ruA * u_int + p_int * A;
    f[2] = rEA * u_int + p_int * u_int * A;
    return f;
  }
  else if (S_R > 0)
  {
    const double factor = A / (S_R - u_int);
    const double rA = (S_R - u_R) * r_R * factor;
    const double ruA = ((S_R - u_R) * r_R * u_R + (p_int - p_R)) * factor;
    const double rEA = ((S_R - u_R) * r_R * E_R - p_R * u_R + p_int * u_int) * factor;

    std::vector<double> f(3, 0.0);
    f[0] = rA * u_int;
    f[1] = ruA * u_int + p_int * A;
    f[2] = rEA * u_int + p_int * u_int * A;
    return f;
  }
  else
  {
    return Euler1Phase::computeFlux(U_R, A, eos_R);
  }
}

void FluxEuler2PhaseDEM::apply1PhaseFlux(
  const std::vector<double> & f_1phase,
  const double & x_times_S,
  const unsigned int & k,
  std::vector<double> & f_2phase) const
{
  f_2phase[3*k + 1] += x_times_S * f_1phase[0];
  f_2phase[3*k + 2] += x_times_S * f_1phase[1];
  f_2phase[3*k + 3] += x_times_S * f_1phase[2];
}

void FluxEuler2PhaseDEM::apply2PhaseFlux(
  const std::vector<double> & f_2phase_add,
  const double & factor,
  std::vector<double> & f_2phase) const
{
  for (unsigned int i = 0; i < f_2phase_add.size(); ++i)
    f_2phase[i] += factor * f_2phase_add[i];
}
