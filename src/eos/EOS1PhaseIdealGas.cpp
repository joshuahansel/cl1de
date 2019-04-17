#include "EOS1PhaseIdealGas.h"
#include "PhysicsConstants.h"

#include <cmath>

EOS1PhaseIdealGas::EOS1PhaseIdealGas(const double & gamma, const double & M)
  : EOS1Phase(),
    _gamma(gamma),
    _M(M),
    _R_sp(PhysicsConstants::gas_constant / _M),
    _cp(_gamma * _R_sp / (_gamma - 1.0)),
    _cv(_cp / _gamma)
{
}

double EOS1PhaseIdealGas::e_from_p_r(const double & p, const double & r) const
{
  return p / ((_gamma - 1) * r);
}

double EOS1PhaseIdealGas::p_from_r_e(const double & r, const double & e) const
{
  return (_gamma - 1) * r * e;
}

double EOS1PhaseIdealGas::T_from_r_e(const double & /*r*/, const double & e) const
{
  return e / _cv;
}

double EOS1PhaseIdealGas::c_from_r_e(const double & r, const double & e) const
{
  const double p = p_from_r_e(r, e);
  return std::sqrt(_gamma * p / r);
}

double EOS1PhaseIdealGas::r_from_p_T(const double & p, const double & T) const
{
  return p / ((_gamma - 1.0) * _cv * T);
}

double EOS1PhaseIdealGas::h_from_p_T(const double & p, const double & T) const
{
  const double r = r_from_p_T(p, T);
  const double e = e_from_p_r(p, r);
  return e + p / r;
}

double EOS1PhaseIdealGas::s_from_p_T(const double & p, const double & T) const
{
  const double n = std::pow(T, _gamma) / std::pow(p, _gamma - 1.0);
  return _cv * std::log(n);
}

double EOS1PhaseIdealGas::e_from_r_h(const double & /*r*/, const double & h) const
{
  return h / _gamma;
}

double EOS1PhaseIdealGas::r_from_h_s(const double & h, const double & s) const
{
  const double p = std::pow(h / (_gamma * _cv), _gamma / (_gamma - 1.0)) *
    std::exp(-s / ((_gamma - 1.0) * _cv));
  const double aux = (s + _cv * std::log(std::pow(p, _gamma - 1.0))) / _cv;
  const double T = std::pow(std::exp(aux), 1.0 / _gamma);
  return r_from_p_T(p, T);
}
