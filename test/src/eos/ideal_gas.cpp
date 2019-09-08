#include <gtest/gtest.h>
#include "EOS1PhaseIdealGas.h"
#include "eos_1phase.h"

TEST(EOS1PhaseIdealGas, consistency)
{
  EOS1PhaseIdealGas eos(1.4, 5.0);
  testEOS1PhaseConsistency(eos);
}
