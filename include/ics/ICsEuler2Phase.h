#ifndef ICsEuler2Phase_H
#define ICsEuler2Phase_H

class EOS1Phase;

class ICsEuler2Phase
{
public:
  ICsEuler2Phase(const EOS1Phase & eos_liq, const EOS1Phase & eos_vap);

  virtual void computeICs(
    const double & x,
    const double & A,
    double & aA_liq,
    double & arA_liq,
    double & aruA_liq,
    double & arEA_liq,
    double & arA_vap,
    double & aruA_vap,
    double & arEA_vap) const = 0;

protected:
  const EOS1Phase & _eos_liq;
  const EOS1Phase & _eos_vap;
};

#endif /* ICsEuler2Phase_H */
