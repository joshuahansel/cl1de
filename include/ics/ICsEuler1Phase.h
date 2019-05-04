#ifndef ICsEuler1Phase_H
#define ICsEuler1Phase_H

class EOS1Phase;

class ICsEuler1Phase
{
public:
  ICsEuler1Phase(const EOS1Phase & eos);

  virtual void computeICs(const double & x, const double & A, double & rA, double & ruA, double & rEA) const = 0;

protected:
  const EOS1Phase & _eos;
};

#endif /* ICsEuler1Phase_H */
