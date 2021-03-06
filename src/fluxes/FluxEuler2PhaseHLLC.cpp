#include "FluxEuler2PhaseHLLC.h"
#include "UtilsEuler2Phase.h"
#include "UtilsNumerics.h"
#include "EOS1Phase.h"
#include "utils.h"

#include <cmath>
#include <sstream>

FluxEuler2PhaseHLLC::FluxEuler2PhaseHLLC(const EOS1Phase & eos_liq, const EOS1Phase & eos_vap)
  : FluxEuler2Phase(eos_liq, eos_vap), _last_region_indices(_n_phases), _a(_n_phases, std::vector<double>(_n_sides))
{
}

void FluxEuler2PhaseHLLC::computeFlux(
  const std::vector<double> & U_L,
  const std::vector<double> & U_R,
  const double & A,
  std::vector<double> & f_L,
  std::vector<double> & f_R) const
{
  const auto f_cons = computeConservativeFlux(U_L, U_R, A);
  const auto f_noncons = Euler2Phase::computeNonConservativePhaseFlux(_u_int, _p_int, A);

  f_L[0] = f_cons[liq][0] + _a[liq][L] * f_noncons[0];
  f_L[1] = f_cons[liq][1] + _a[liq][L] * f_noncons[1];
  f_L[2] = f_cons[liq][2] + _a[liq][L] * f_noncons[2];
  f_L[3] = f_cons[liq][3] + _a[liq][L] * f_noncons[3];
  f_L[4] = f_cons[vap][1] + _a[vap][L] * f_noncons[1];
  f_L[5] = f_cons[vap][2] + _a[vap][L] * f_noncons[2];
  f_L[6] = f_cons[vap][3] + _a[vap][L] * f_noncons[3];

  f_R[0] = f_cons[liq][0] + _a[liq][R] * f_noncons[0];
  f_R[1] = f_cons[liq][1] + _a[liq][R] * f_noncons[1];
  f_R[2] = f_cons[liq][2] + _a[liq][R] * f_noncons[2];
  f_R[3] = f_cons[liq][3] + _a[liq][R] * f_noncons[3];
  f_R[4] = f_cons[vap][1] + _a[vap][R] * f_noncons[1];
  f_R[5] = f_cons[vap][2] + _a[vap][R] * f_noncons[2];
  f_R[6] = f_cons[vap][3] + _a[vap][R] * f_noncons[3];
}

std::vector<std::vector<double>> FluxEuler2PhaseHLLC::computeConservativeFlux(
  const std::vector<double> & U_L,
  const std::vector<double> & U_R,
  const double & A) const
{
  const std::vector<std::vector<double>> U {U_L, U_R};

  // unpack data
  std::vector<std::vector<double>> aA(_n_phases, std::vector<double>(_n_sides));
  std::vector<std::vector<double>> arA(_n_phases, std::vector<double>(_n_sides));
  std::vector<std::vector<double>> aruA(_n_phases, std::vector<double>(_n_sides));
  std::vector<std::vector<double>> arEA(_n_phases, std::vector<double>(_n_sides));
  for (unsigned int i = 0; i < _n_sides; i++)
  {
    aA[liq][i] = U[i][0];
    arA[liq][i] = U[i][1];
    aruA[liq][i] = U[i][2];
    arEA[liq][i] = U[i][3];
    arA[vap][i] = U[i][4];
    aruA[vap][i] = U[i][5];
    arEA[vap][i] = U[i][6];

    aA[vap][i] = (1.0 - aA[liq][i] / A) * A;
  }

  // compute and store primitive variables
  std::vector<std::vector<double>> r(_n_phases, std::vector<double>(_n_sides));
  std::vector<std::vector<double>> u(_n_phases, std::vector<double>(_n_sides));
  std::vector<std::vector<double>> E(_n_phases, std::vector<double>(_n_sides));
  std::vector<std::vector<double>> e(_n_phases, std::vector<double>(_n_sides));
  std::vector<std::vector<double>> p(_n_phases, std::vector<double>(_n_sides));
  std::vector<std::vector<double>> c(_n_phases, std::vector<double>(_n_sides));
  std::vector<std::vector<std::vector<double>>> W(
    _n_phases, std::vector<std::vector<double>>(_n_sides, std::vector<double>(Euler2Phase::n_local_prim_var)));
  for (unsigned int i = 0; i < _n_sides; i++)
  {
    for (unsigned int k = 0; k < _n_phases; k++)
    {
      _a[k][i] = aA[k][i] / A;
      r[k][i] = arA[k][i] / (_a[k][i] * A);
      u[k][i] = aruA[k][i] / arA[k][i];
      E[k][i] = arEA[k][i] / arA[k][i];
      e[k][i] = E[k][i] - 0.5 * u[k][i] * u[k][i];
      p[k][i] = _eos[k]->p_from_r_e(r[k][i], e[k][i]);
      c[k][i] = _eos[k]->c_from_r_e(r[k][i], e[k][i]);

      W[k][i][Euler2Phase::ia] = _a[k][i];
      W[k][i][Euler2Phase::ir] = r[k][i];
      W[k][i][Euler2Phase::iu] = u[k][i];
      W[k][i][Euler2Phase::ip] = p[k][i];
      W[k][i][Euler2Phase::iE] = E[k][i];
    }
  }

  // compute tentative wave speeds
  std::vector<std::vector<double>> St(_n_phases, std::vector<double>(_n_sides));
  for (unsigned int k = 0; k < _n_phases; k++)
  {
    St[k][L] = std::min(u[k][L] - c[k][L], u[k][R] - c[k][R]);
    St[k][R] = std::max(u[k][L] + c[k][L], u[k][R] + c[k][R]);
  }

  // compute interfacial velocity and pressure
  const unsigned int kk = _a[liq][L] < _a[liq][R] ? vap : liq;
  const unsigned int jj = _a[liq][L] < _a[liq][R] ? liq : vap;
  _u_int = ((r[jj][R] * u[jj][R] * (u[jj][R] - St[jj][R]) + p[jj][R]) - (r[kk][L] * u[kk][L] * (u[kk][L] - St[kk][L]) + p[kk][L]))
    / ((r[jj][R] * (u[jj][R] - St[jj][R])) - (r[kk][L] * (u[kk][L] - St[kk][L])));
  _p_int = r[jj][R] * (u[jj][R] - St[jj][R]) * (u[jj][R] - _u_int) + p[jj][R];

  // compute limited wave speeds
  std::vector<std::vector<double>> S(_n_phases, std::vector<double>(_n_sides));
  std::vector<double> SM(_n_phases);
  for (unsigned int k = 0; k < _n_phases; k++)
  {
    S[k][L] = _u_int < St[k][L] + 0.1 * std::abs(St[k][L]) ? std::min(St[liq][L], St[vap][L]) : St[k][L];
    S[k][R] = _u_int > St[k][R] - 0.1 * std::abs(St[k][R]) ? std::max(St[liq][R], St[vap][R]) : St[k][R];
    SM[k] = (
        (_a[k][R] * r[k][R] * u[k][R] * (u[k][R] - S[k][R]) + _a[k][R] * p[k][R])
      - (_a[k][L] * r[k][L] * u[k][L] * (u[k][L] - S[k][L]) + _a[k][L] * p[k][L])
      + (_a[k][L] - _a[k][R]) * _p_int
    ) / (
        (_a[k][R] * r[k][R] * (u[k][R] - S[k][R]))
      - (_a[k][L] * r[k][L] * (u[k][L] - S[k][L]))
    );
  }

  // compute conservative flux
  std::vector<std::vector<double>> f_cons(_n_phases, std::vector<double>(Euler2Phase::n_local_eq));
  for (unsigned int k = 0; k < _n_phases; k++)
  {
    if (_u_int < SM[k])
    {
      if (S[k][L] > 0)
      {
        f_cons[k] = Euler2Phase::computeConservativePhaseFlux(W[k][L], _u_int, _p_int, A);

        _last_region_indices[k] = 0;
      }
      else if (_u_int > 0)
      {
        const auto UL = solutionSupersonic(W[k][L], A);
        const auto Us = solutionSubsonic(W[k][L], S[k][L], SM[k], A);

        f_cons[k] = Euler2Phase::computeConservativePhaseFlux(W[k][L], _u_int, _p_int, A);
        applyRankineHugoniot(f_cons[k], S[k][L], UL, Us);

        _last_region_indices[k] = 1;
      }
      else if (SM[k] > 0)
      {
        const auto UR = solutionSupersonic(W[k][R], A);
        const auto Us = solutionSubsonic(W[k][R], S[k][R], SM[k], A);
        const auto Uss = solutionSubsonicInterfacialLeft(W[k][L], W[k][R], S[k][L], S[k][R], SM[k], _p_int, A);

        f_cons[k] = Euler2Phase::computeConservativePhaseFlux(W[k][R], _u_int, _p_int, A);
        applyRankineHugoniot(f_cons[k], S[k][R], UR, Us);
        applyRankineHugoniot(f_cons[k], SM[k], Us, Uss);

        _last_region_indices[k] = 2;
      }
      else if (S[k][R] > 0)
      {
        const auto UR = solutionSupersonic(W[k][R], A);
        const auto Us = solutionSubsonic(W[k][R], S[k][R], SM[k], A);

        f_cons[k] = Euler2Phase::computeConservativePhaseFlux(W[k][R], _u_int, _p_int, A);
        applyRankineHugoniot(f_cons[k], S[k][R], UR, Us);

        _last_region_indices[k] = 3;
      }
      else if (S[k][R] <= 0)
      {
        f_cons[k] = Euler2Phase::computeConservativePhaseFlux(W[k][R], _u_int, _p_int, A);

        _last_region_indices[k] = 4;
      }
      else
      {
        std::stringstream ss;
        ss << "Riemann solver encountered a NaN:\n";
        ss << "  k = " << k << "\n";
        ss << "  u_int = " << _u_int << "\n";
        ss << "  SL = " << S[k][L] << "\n";
        ss << "  SM = " << SM[k] << "\n";
        ss << "  SR = " << S[k][R] << "\n";
        for (unsigned int i = 0; i < _n_sides; i++)
        {
          ss << "  i = " << i << ":\n";
          ss << "    a = " << _a[k][i] << "\n";
          ss << "    r = " << r[k][i] << "\n";
          ss << "    u = " << u[k][i] << "\n";
          ss << "    E = " << E[k][i] << "\n";
          ss << "    e = " << e[k][i] << "\n";
          ss << "    p = " << p[k][i] << "\n";
          ss << "    c = " << c[k][i] << "\n";
        }
        throwError(ss.str());
      }
    }
    else if (_u_int >= SM[k])
    {
      if (S[k][L] > 0)
      {
        f_cons[k] = Euler2Phase::computeConservativePhaseFlux(W[k][L], _u_int, _p_int, A);

        _last_region_indices[k] = 5;
      }
      else if (SM[k] > 0)
      {
        const auto UL = solutionSupersonic(W[k][L], A);
        const auto Us = solutionSubsonic(W[k][L], S[k][L], SM[k], A);

        f_cons[k] = Euler2Phase::computeConservativePhaseFlux(W[k][L], _u_int, _p_int, A);
        applyRankineHugoniot(f_cons[k], S[k][L], UL, Us);

        _last_region_indices[k] = 6;
      }
      else if (_u_int > 0)
      {
        const auto UL = solutionSupersonic(W[k][L], A);
        const auto Us = solutionSubsonic(W[k][L], S[k][L], SM[k], A);
        const auto Uss = solutionSubsonicInterfacialLeft(W[k][R], W[k][L], S[k][R], S[k][L], SM[k], _p_int, A);

        f_cons[k] = Euler2Phase::computeConservativePhaseFlux(W[k][L], _u_int, _p_int, A);
        applyRankineHugoniot(f_cons[k], S[k][L], UL, Us);
        applyRankineHugoniot(f_cons[k], SM[k], Us, Uss);

        _last_region_indices[k] = 7;
      }
      else if (S[k][R] > 0)
      {
        const auto UR = solutionSupersonic(W[k][R], A);
        const auto Us = solutionSubsonic(W[k][R], S[k][R], SM[k], A);

        f_cons[k] = Euler2Phase::computeConservativePhaseFlux(W[k][R], _u_int, _p_int, A);
        applyRankineHugoniot(f_cons[k], S[k][R], UR, Us);

        _last_region_indices[k] = 8;
      }
      else if (S[k][R] <= 0)
      {
        f_cons[k] = Euler2Phase::computeConservativePhaseFlux(W[k][R], _u_int, _p_int, A);

        _last_region_indices[k] = 9;
      }
      else
      {
        std::stringstream ss;
        ss << "Riemann solver encountered a NaN:\n";
        ss << "  k = " << k << "\n";
        ss << "  u_int = " << _u_int << "\n";
        ss << "  SL = " << S[k][L] << "\n";
        ss << "  SM = " << SM[k] << "\n";
        ss << "  SR = " << S[k][R] << "\n";
        for (unsigned int i = 0; i < _n_sides; i++)
        {
          ss << "  i = " << i << ":\n";
          ss << "    a = " << _a[k][i] << "\n";
          ss << "    r = " << r[k][i] << "\n";
          ss << "    u = " << u[k][i] << "\n";
          ss << "    E = " << E[k][i] << "\n";
          ss << "    e = " << e[k][i] << "\n";
          ss << "    p = " << p[k][i] << "\n";
          ss << "    c = " << c[k][i] << "\n";
        }
        throwError(ss.str());
      }
    }
    else
    {
      std::stringstream ss;
      ss << "Riemann solver encountered a NaN:\n";
      ss << "  k = " << k << "\n";
      ss << "  u_int = " << _u_int << "\n";
      ss << "  SL = " << S[k][L] << "\n";
      ss << "  SM = " << SM[k] << "\n";
      ss << "  SR = " << S[k][R] << "\n";
      for (unsigned int i = 0; i < _n_sides; i++)
      {
        ss << "  i = " << i << ":\n";
        ss << "    a = " << _a[k][i] << "\n";
        ss << "    r = " << r[k][i] << "\n";
        ss << "    u = " << u[k][i] << "\n";
        ss << "    E = " << E[k][i] << "\n";
        ss << "    e = " << e[k][i] << "\n";
        ss << "    p = " << p[k][i] << "\n";
        ss << "    c = " << c[k][i] << "\n";
      }
      throwError(ss.str());
    }
  }

  return f_cons;
}

std::vector<double> FluxEuler2PhaseHLLC::solutionSupersonic(
  const std::vector<double> & W, const double & A) const
{
  const double & a = W[Euler2Phase::ia];
  const double & r = W[Euler2Phase::ir];
  const double & u = W[Euler2Phase::iu];
  const double & E = W[Euler2Phase::iE];

  std::vector<double> U(Euler2Phase::n_local_eq);
  U[0] = a * A;
  U[1] = a * r * A;
  U[2] = a * r * u * A;
  U[3] = a * r * E * A;

  return U;
}

std::vector<double> FluxEuler2PhaseHLLC::solutionSubsonic(
  const std::vector<double> & W,
  const double & S,
  const double & SM,
  const double & A
) const
{
  const double & a = W[Euler2Phase::ia];
  const double & r = W[Euler2Phase::ir];
  const double & u = W[Euler2Phase::iu];
  const double & p = W[Euler2Phase::ip];
  const double & E = W[Euler2Phase::iE];

  const double as = a;
  const double rs = r * (u - S) / (SM - S);
  const double us = SM;
  const double ps = p + r * (u - S) * (u - SM);
  const double Es = E + (p * u - ps * SM) / (r * (u - S));

  std::vector<double> Us(Euler2Phase::n_local_eq);
  Us[0] = as * A;
  Us[1] = as * rs * A;
  Us[2] = as * rs * us * A;
  Us[3] = as * rs * Es * A;

  return Us;
}

std::vector<double> FluxEuler2PhaseHLLC::solutionSubsonicInterfacialLeft(
  const std::vector<double> & WL,
  const std::vector<double> & WR,
  const double & SL,
  const double & SR,
  const double & SM,
  const double & p_int,
  const double & A
) const
{
  const double & aL = WL[Euler2Phase::ia];
  const double & rL = WL[Euler2Phase::ir];
  const double & uL = WL[Euler2Phase::iu];
  const double & pL = WL[Euler2Phase::ip];
  const double & EL = WL[Euler2Phase::iE];

  const double asL = aL;
  const double asR = WR[Euler2Phase::ia];
  const double rsL = rL * (uL - SL) / (SM - SL);
  const double psL = pL + rL * (uL - SL) * (uL - SM);
  const double EsL = EL + (pL * uL - psL * SM) / (rL * (uL - SL));

  const double ass = asR;
  const double rss = (asL * rsL) / ass;
  const double uss = SM;
  const double Ess = EsL - (asR - asL) / (asL * rsL) * p_int;

  std::vector<double> Uss(Euler2Phase::n_local_eq);
  Uss[0] = ass * A;
  Uss[1] = ass * rss * A;
  Uss[2] = ass * rss * uss * A;
  Uss[3] = ass * rss * Ess * A;

  return Uss;
}

void FluxEuler2PhaseHLLC::applyRankineHugoniot(
  std::vector<double> & flux,
  const double & S,
  const std::vector<double> & U1,
  const std::vector<double> & U2
) const
{
  for (unsigned int i = 0; i < Euler2Phase::n_local_eq; i++)
    flux[i] += S * (U2[i] - U1[i]);
}
