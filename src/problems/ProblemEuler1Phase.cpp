#include "ProblemEuler1Phase.h"
#include "Function.h"
#include "EOS1Phase.h"
#include "BCEuler1Phase.h"

ProblemEuler1Phase::ProblemEuler1Phase() : Problem()
{
}

void ProblemEuler1Phase::setAreaFunction(const Function & function)
{
  _A_fn = &function;
}

void ProblemEuler1Phase::setICDensity(const Function & function)
{
  _r_ic_fn = &function;
}

void ProblemEuler1Phase::setICVelocity(const Function & function)
{
  _u_ic_fn = &function;
}

void ProblemEuler1Phase::setICPressure(const Function & function)
{
  _p_ic_fn = &function;
}

void ProblemEuler1Phase::setEOS(const EOS1Phase & eos)
{
  _eos = &eos;
}

void ProblemEuler1Phase::setBCLeft(const BCEuler1Phase & bc)
{
  _bc_left = &bc;
}

void ProblemEuler1Phase::setBCRight(const BCEuler1Phase & bc)
{
  _bc_right = &bc;
}
