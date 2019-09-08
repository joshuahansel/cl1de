#include "eos_1phase.h"
#include "TestingUtils.h"

void testEOS1PhaseConsistency(const EOS1Phase & eos)
{
  const double p = 1e5;
  const double T = 300;

  const double r_from_p_T = eos.r_from_p_T(p, T);
  const double h_from_p_T = eos.h_from_p_T(p, T);
  const double s_from_p_T = eos.s_from_p_T(p, T);

  const double r_from_h_s = eos.r_from_h_s(h_from_p_T, s_from_p_T);
  ABS_TEST(r_from_p_T, r_from_h_s, 1e-12);

  const double e_from_p_r = eos.e_from_p_r(p, r_from_p_T);
  const double e_from_r_h = eos.e_from_r_h(r_from_p_T, h_from_p_T);
  ABS_TEST(e_from_p_r, e_from_r_h, 1e-11);

  const double p_from_r_e = eos.p_from_r_e(r_from_p_T, e_from_p_r);
  ABS_TEST(p_from_r_e, p, 1e-12);

  const double T_from_r_e = eos.T_from_r_e(r_from_p_T, e_from_p_r);
  ABS_TEST(T_from_r_e, T, 1e-12);
}
