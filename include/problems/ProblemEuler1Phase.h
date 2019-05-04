#ifndef ProblemEuler1Phase_H
#define ProblemEuler1Phase_H

#include "Problem.h"

class Function;
class ICsEuler1Phase;
class EOS1Phase;
class BCEuler1Phase;

class ProblemEuler1Phase : public Problem
{
public:
  ProblemEuler1Phase();

  void setAreaFunction(const Function & function);
  void setICs(const ICsEuler1Phase & ics);
  void setEOS(const EOS1Phase & eos);
  void setBCLeft(const BCEuler1Phase & bc);
  void setBCRight(const BCEuler1Phase & bc);

  const Function & getAreaFunction() const {return *_A_fn;}
  const ICsEuler1Phase & getICs() const {return *_ics;}
  const EOS1Phase & getEOS() const {return *_eos;}
  const BCEuler1Phase & getBCLeft() const {return *_bc_left;}
  const BCEuler1Phase & getBCRight() const {return *_bc_right;}

protected:
  const Function * _A_fn;
  const ICsEuler1Phase * _ics;
  const EOS1Phase * _eos;
  const BCEuler1Phase * _bc_left;
  const BCEuler1Phase * _bc_right;
};

#endif /* ProblemEuler1Phase_H */
