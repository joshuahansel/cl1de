#ifndef EOS1PhaseIdealGas_H
#define EOS1PhaseIdealGas_H

#include "EOS1Phase.h"

class EOS1PhaseIdealGas : public EOS1Phase
{
public:
  EOS1PhaseIdealGas(const double & gamma, const double & M);

  virtual double e_from_p_r(const double & p, const double & r) const override;

  virtual double p_from_r_e(const double & r, const double & e) const override;
  virtual double T_from_r_e(const double & r, const double & e) const override;
  virtual double c_from_r_e(const double & r, const double & e) const override;

  virtual double r_from_p_T(const double & p, const double & T) const override;
  virtual double h_from_p_T(const double & p, const double & T) const override;
  virtual double s_from_p_T(const double & p, const double & T) const override;

  virtual double e_from_r_h(const double & r, const double & h) const override;

  virtual double r_from_h_s(const double & h, const double & s) const override;

protected:
  const double _gamma;
  const double _M;
  const double _R_sp;
  const double _cp;
  const double _cv;
};

#endif /* EOS1PhaseIdealGas_H */
