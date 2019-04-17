#ifndef EOS1Phase_H
#define EOS1Phase_H

class EOS1Phase
{
public:
  EOS1Phase();

  virtual double e_from_p_r(const double & p, const double & r) const = 0;

  virtual double p_from_r_e(const double & r, const double & e) const = 0;
  virtual double T_from_r_e(const double & r, const double & e) const = 0;
  virtual double c_from_r_e(const double & r, const double & e) const = 0;

  virtual double r_from_p_T(const double & p, const double & T) const = 0;
  virtual double h_from_p_T(const double & p, const double & T) const = 0;
  virtual double s_from_p_T(const double & p, const double & T) const = 0;

  virtual double e_from_r_h(const double & r, const double & h) const = 0;

  virtual double r_from_h_s(const double & h, const double & s) const = 0;
};

#endif /* EOS1Phase_H */
