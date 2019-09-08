#include "ProblemEuler1Phase.h"
#include "Function.h"
#include "ICsEuler1Phase.h"
#include "EOS1Phase.h"
#include "BCEuler1Phase.h"

ProblemEuler1Phase::ProblemEuler1Phase() : Problem(), _gravity(0.0)
{
}

void ProblemEuler1Phase::setAreaFunction(const Function & function)
{
  _A_fn = &function;
}

void ProblemEuler1Phase::setICs(const ICsEuler1Phase & ics)
{
  _ics = &ics;
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

void ProblemEuler1Phase::setGravity(const double & gravity)
{
  _gravity = gravity;
}
