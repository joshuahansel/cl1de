#include "ProblemEuler2Phase.h"
#include "Function.h"
#include "ICsEuler2Phase.h"
#include "EOS1Phase.h"
#include "BCEuler2Phase.h"

ProblemEuler2Phase::ProblemEuler2Phase() : Problem(), _gravity(0.0), _a_int(0.0)
{
}

void ProblemEuler2Phase::setAreaFunction(const Function & function)
{
  _A_fn = &function;
}

void ProblemEuler2Phase::setICs(const ICsEuler2Phase & ics)
{
  _ics = &ics;
}

void ProblemEuler2Phase::setEOS(const EOS1Phase & eos_liq, const EOS1Phase & eos_vap)
{
  _eos_liq = &eos_liq;
  _eos_vap = &eos_vap;
}

void ProblemEuler2Phase::setBCLeft(const BCEuler2Phase & bc)
{
  _bc_left = &bc;
}

void ProblemEuler2Phase::setBCRight(const BCEuler2Phase & bc)
{
  _bc_right = &bc;
}

void ProblemEuler2Phase::setGravity(const double & gravity)
{
  _gravity = gravity;
}

void ProblemEuler2Phase::setInterfacialAreaDensity(const double & a_int)
{
  _a_int = a_int;
}
