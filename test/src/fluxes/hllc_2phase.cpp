#include "FluxEuler2PhaseHLLC.h"
#include "EOS1PhaseIdealGas.h"
#include "PhysicsConstants.h"
#include "TestingUtils.h"

#include <gtest/gtest.h>

namespace TestFluxEuler2PhaseHLLC
{
std::vector<double>
computeConservativeSolution(
  const std::vector<double> & W, const double & A, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap)
{
  const double & a_liq = W[0];
  const double & p_liq = W[1];
  const double & T_liq = W[2];
  const double & u_liq = W[3];
  const double & p_vap = W[4];
  const double & T_vap = W[5];
  const double & u_vap = W[6];

  const double r_liq = eos_liq.r_from_p_T(p_liq, T_liq);
  const double e_liq = eos_liq.e_from_p_r(p_liq, r_liq);
  const double E_liq = e_liq + 0.5 * u_liq * u_liq;

  const double a_vap = 1.0 - a_liq;
  const double r_vap = eos_vap.r_from_p_T(p_vap, T_vap);
  const double e_vap = eos_vap.e_from_p_r(p_vap, r_vap);
  const double E_vap = e_vap + 0.5 * u_vap * u_vap;

  std::vector<double> U(7, 0.0);
  U[0] = a_liq * A;
  U[1] = a_liq * r_liq * A;
  U[2] = a_liq * r_liq * u_liq * A;
  U[3] = a_liq * r_liq * E_liq * A;
  U[4] = a_vap * r_vap * A;
  U[5] = a_vap * r_vap * u_vap * A;
  U[6] = a_vap * r_vap * E_vap * A;

  return U;
}
}

TEST(FluxEuler2PhaseHLLC, symmetry)
{
  const double gamma_liq = 2.35;
  const double gamma_vap = 1.43;
  const double R_liq = 2451.6;
  const double R_vap = 447.2;
  const double M_liq = PhysicsConstants::gas_constant / R_liq;
  const double M_vap = PhysicsConstants::gas_constant / R_vap;

  EOS1PhaseIdealGas eos_liq(gamma_liq, M_liq);
  EOS1PhaseIdealGas eos_vap(gamma_vap, M_vap);

  FluxEuler2PhaseHLLC flux(eos_liq, eos_vap);

  const double A = 0.3;

  std::vector<std::pair<std::vector<double>, std::vector<double>>> W_pairs;
  {
    const std::vector<double> W_L{0.3, 1e5, 300, 1.5, 1.2e5, 305, 1.2};
    const std::vector<double> W_R{0.6, 2e5, 310, 1.3, 2.1e5, 280, 0.7};
    W_pairs.push_back(std::pair<std::vector<double>, std::vector<double>>(W_L, W_R));
  }
  {
    const std::vector<double> W_L{0.3, 1e5, 300, 1500, 1.2e5, 305, -1000};
    const std::vector<double> W_R{0.6, 2e5, 310, 1400, 2.1e5, 280, 0.7};
    W_pairs.push_back(std::pair<std::vector<double>, std::vector<double>>(W_L, W_R));
  }
  {
    const std::vector<double> W_L{0.3, 5e5, 300, -1500, 1.2e5, 305, -1450};
    const std::vector<double> W_R{0.6, 2e5, 310, -1400, 2.1e5, 280, -1600};
    W_pairs.push_back(std::pair<std::vector<double>, std::vector<double>>(W_L, W_R));
  }

  std::vector<double> f_L(7, 0.0), f_R(7, 0.0), f_L_mirror(7, 0.0), f_R_mirror(7, 0.0);
  std::set<unsigned int> flux_regions;
  for (const auto & W_pair : W_pairs)
  {
    const auto & W_L = W_pair.first;
    const auto & W_R = W_pair.second;

    const auto U_L = TestFluxEuler2PhaseHLLC::computeConservativeSolution(W_L, A, eos_liq, eos_vap);
    const auto U_R = TestFluxEuler2PhaseHLLC::computeConservativeSolution(W_R, A, eos_liq, eos_vap);

    auto U_L_mirror = U_L;
    U_L_mirror[2] *= -1.0;
    U_L_mirror[5] *= -1.0;

    auto U_R_mirror = U_R;
    U_R_mirror[2] *= -1.0;
    U_R_mirror[5] *= -1.0;

    flux.computeFlux(U_L, U_R, A, f_L, f_R);
    auto indices = flux.getLastRegionIndices();
    std::copy(indices.begin(), indices.end(), std::inserter(flux_regions, flux_regions.end()));

    flux.computeFlux(U_R_mirror, U_L_mirror, A, f_R_mirror, f_L_mirror);
    indices = flux.getLastRegionIndices();
    std::copy(indices.begin(), indices.end(), std::inserter(flux_regions, flux_regions.end()));

    REL_TEST(f_L_mirror[0], -f_L[0], 1e-12);
    REL_TEST(f_L_mirror[1], -f_L[1], 1e-12);
    REL_TEST(f_L_mirror[2],  f_L[2], 1e-12);
    REL_TEST(f_L_mirror[3], -f_L[3], 1e-12);
    REL_TEST(f_L_mirror[4], -f_L[4], 1e-12);
    REL_TEST(f_L_mirror[5],  f_L[5], 1e-12);
    REL_TEST(f_L_mirror[6], -f_L[6], 1e-12);

    REL_TEST(f_R_mirror[0], -f_R[0], 1e-12);
    REL_TEST(f_R_mirror[1], -f_R[1], 1e-12);
    REL_TEST(f_R_mirror[2],  f_R[2], 1e-12);
    REL_TEST(f_R_mirror[3], -f_R[3], 1e-12);
    REL_TEST(f_R_mirror[4], -f_R[4], 1e-12);
    REL_TEST(f_R_mirror[5],  f_R[5], 1e-12);
    REL_TEST(f_R_mirror[6], -f_R[6], 1e-12);
  }

  EXPECT_TRUE(flux_regions.size() == 10);
}
