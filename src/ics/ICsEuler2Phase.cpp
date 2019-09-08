#include "ICsEuler2Phase.h"

ICsEuler2Phase::ICsEuler2Phase(const EOS1Phase & eos_liq, const EOS1Phase & eos_vap)
  : _eos_liq(eos_liq), _eos_vap(eos_vap)
{
}
