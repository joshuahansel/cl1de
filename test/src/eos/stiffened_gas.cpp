#include <gtest/gtest.h>
#include "TestingUtils.h"
#include "EOS1PhaseStiffenedGas.h"
#include "EOS1PhaseIdealGas.h"
#include "PhysicsConstants.h"
#include "eos_1phase.h"

TEST(EOS1PhaseStiffenedGas, consistency)
{
  EOS1PhaseStiffenedGas eos(1.4, 50.0, 1e3, 5.0, 8.0);
  testEOS1PhaseConsistency(eos);
}

TEST(EOS1PhaseStiffenedGas, ideal_gas_equivalency)
{
  const double gamma = 1.4;
  const double M = 0.001;
  const double R_sp = PhysicsConstants::gas_constant / M;
  const double cp = gamma * R_sp / (gamma - 1.0);
  const double cv = cp / gamma;

  EOS1PhaseStiffenedGas eos_stiffened(gamma, cv, 0, 0, 0);
  EOS1PhaseIdealGas eos_ideal(gamma, M);

  const double p = 1e5;
  const double T = 300;
  const double r = eos_ideal.r_from_p_T(p, T);
  const double h = eos_ideal.h_from_p_T(p, T);
  const double s = eos_ideal.s_from_p_T(p, T);
  const double e = eos_ideal.e_from_p_r(p, r);
  const double c = eos_ideal.c_from_r_e(r, e);

  ABS_TEST(eos_stiffened.r_from_p_T(p, T), r, 1e-12);
  ABS_TEST(eos_stiffened.h_from_p_T(p, T), h, 1e-12);
  ABS_TEST(eos_stiffened.s_from_p_T(p, T), s, 1e-12);
  ABS_TEST(eos_stiffened.e_from_p_r(p, r), e, 1e-12);
  ABS_TEST(eos_stiffened.c_from_r_e(r, e), c, 1e-12);
}
