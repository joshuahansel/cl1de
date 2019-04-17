#ifndef ProblemEuler1Phase_H
#define ProblemEuler1Phase_H

#include "Problem.h"

class Function;
class EOS1Phase;
class BCEuler1Phase;

class ProblemEuler1Phase : public Problem
{
public:
  ProblemEuler1Phase();

  void setAreaFunction(const Function & function);
  void setICDensity(const Function & function);
  void setICVelocity(const Function & function);
  void setICPressure(const Function & function);
  void setEOS(const EOS1Phase & eos);
  void setBCLeft(const BCEuler1Phase & bc);
  void setBCRight(const BCEuler1Phase & bc);

  const Function & getAreaFunction() const {return *_A_fn;}
  const Function & getICDensity() const {return *_r_ic_fn;}
  const Function & getICVelocity() const {return *_u_ic_fn;}
  const Function & getICPressure() const {return *_p_ic_fn;}
  const EOS1Phase & getEOS() const {return *_eos;}
  const BCEuler1Phase & getBCLeft() const {return *_bc_left;}
  const BCEuler1Phase & getBCRight() const {return *_bc_right;}

protected:
  const Function * _A_fn;
  const Function * _r_ic_fn;
  const Function * _u_ic_fn;
  const Function * _p_ic_fn;
  const EOS1Phase * _eos;
  const BCEuler1Phase * _bc_left;
  const BCEuler1Phase * _bc_right;
};

#endif /* ProblemEuler1Phase_H */
