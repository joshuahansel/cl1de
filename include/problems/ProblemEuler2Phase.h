#ifndef ProblemEuler2Phase_H
#define ProblemEuler2Phase_H

#include "Problem.h"

class Function;
class ICsEuler2Phase;
class EOS1Phase;
class BCEuler2Phase;

class ProblemEuler2Phase : public Problem
{
public:
  ProblemEuler2Phase();

  void setAreaFunction(const Function & function);
  void setICs(const ICsEuler2Phase & ics);
  void setEOS(const EOS1Phase & eos_liq, const EOS1Phase & eos_vap);
  void setBCLeft(const BCEuler2Phase & bc);
  void setBCRight(const BCEuler2Phase & bc);
  void setGravity(const double & gravity);
  void setInterfacialAreaDensity(const double & a_int);

  const Function & getAreaFunction() const {return *_A_fn;}
  const ICsEuler2Phase & getICs() const {return *_ics;}
  const EOS1Phase & getEOSLiquid() const {return *_eos_liq;}
  const EOS1Phase & getEOSVapor() const {return *_eos_vap;}
  const BCEuler2Phase & getBCLeft() const {return *_bc_left;}
  const BCEuler2Phase & getBCRight() const {return *_bc_right;}
  const double & getGravity() const {return _gravity;}
  const double & getInterfacialAreaDensity() const {return _a_int;}

protected:
  const Function * _A_fn;
  const ICsEuler2Phase * _ics;
  const EOS1Phase * _eos_liq;
  const EOS1Phase * _eos_vap;
  const BCEuler2Phase * _bc_left;
  const BCEuler2Phase * _bc_right;
  double _gravity;
  double _a_int;
};

#endif /* ProblemEuler2Phase_H */
