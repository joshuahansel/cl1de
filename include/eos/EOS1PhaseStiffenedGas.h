#ifndef EOS1PhaseStiffenedGas_H
#define EOS1PhaseStiffenedGas_H

#include "EOS1Phase.h"

class EOS1PhaseStiffenedGas : public EOS1Phase
{
public:
  EOS1PhaseStiffenedGas(
    const double & gamma,
    const double & cv,
    const double & p_inf,
    const double & q,
    const double & q_prime);

  virtual double e_from_p_r(const double & p, const double & r) const override;

  virtual double p_from_r_e(const double & r, const double & e) const override;
  virtual double T_from_r_e(const double & r, const double & e) const override;
  virtual double c_from_r_e(const double & r, const double & e) const override;

  virtual double r_from_p_T(const double & p, const double & T) const override;
  virtual double h_from_p_T(const double & p, const double & T) const override;
  virtual double s_from_p_T(const double & p, const double & T) const override;

  virtual double e_from_r_h(const double & r, const double & h) const override;

  virtual double r_from_h_s(const double & h, const double & s) const override;

  const double & getGamma() const {return _gamma;}
  const double & getReferencePressure() const {return _p_inf;}

protected:
  const double _gamma;
  const double _cv;
  const double _p_inf;
  const double _q;
  const double _q_prime;
};

#endif /* EOS1PhaseStiffenedGas_H */
