#include "EOS1PhaseStiffenedGas.h"

#include <cmath>

EOS1PhaseStiffenedGas::EOS1PhaseStiffenedGas(
  const double & gamma,
  const double & cv,
  const double & p_inf,
  const double & q,
  const double & q_prime)
  : EOS1Phase(),
    _gamma(gamma),
    _cv(cv),
    _p_inf(p_inf),
    _q(q),
    _q_prime(q_prime)
{
}

double EOS1PhaseStiffenedGas::e_from_p_r(const double & p, const double & r) const
{
  return (p + _gamma * _p_inf) / ((_gamma - 1) * r) + _q;
}

double EOS1PhaseStiffenedGas::p_from_r_e(const double & r, const double & e) const
{
  return (_gamma - 1) * r * (e  - _q) - _gamma * _p_inf;
}

double EOS1PhaseStiffenedGas::T_from_r_e(const double & r, const double & e) const
{
  return (e - _q - _p_inf / r) / _cv;
}

double EOS1PhaseStiffenedGas::c_from_r_e(const double & r, const double & e) const
{
  const double p = p_from_r_e(r, e);
  return std::sqrt(_gamma * (p + _p_inf) / r);
}

double EOS1PhaseStiffenedGas::r_from_p_T(const double & p, const double & T) const
{
  return (p + _p_inf) / ((_gamma - 1.0) * _cv * T);
}

double EOS1PhaseStiffenedGas::h_from_p_T(const double & /*p*/, const double & T) const
{
  return _gamma * _cv * T + _q;
}

double EOS1PhaseStiffenedGas::s_from_p_T(const double & p, const double & T) const
{
  const double n = std::pow(T, _gamma) / std::pow(p + _p_inf, _gamma - 1.0);
  return _cv * std::log(n) + _q_prime;
}

double EOS1PhaseStiffenedGas::e_from_r_h(const double & r, const double & h) const
{
  return (h + (_gamma - 1.0) * _q + _gamma * _p_inf / r) / _gamma;
}

double EOS1PhaseStiffenedGas::r_from_h_s(const double & h, const double & s) const
{
  const double p = std::pow((h - _q) / (_gamma * _cv), _gamma / (_gamma - 1.0)) *
    std::exp((_q_prime - s) / ((_gamma - 1.0) * _cv)) - _p_inf;

  const double aux = (s - _q_prime + _cv * std::log(std::pow(p + _p_inf, _gamma - 1.0))) / _cv;
  const double T = std::pow(std::exp(aux), 1.0 / _gamma);
  return r_from_p_T(p, T);
}
