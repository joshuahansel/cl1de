#include "FluxEuler2Phase.h"
#include "EOS1Phase.h"

FluxEuler2Phase::FluxEuler2Phase(const EOS1Phase & eos_liq, const EOS1Phase & eos_vap)
  : _eos_liq(eos_liq), _eos_vap(eos_vap), _eos({&eos_liq, &eos_vap})
{
}
