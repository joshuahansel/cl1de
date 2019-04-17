#include "UtilsNumerics.h"

#include <cstdlib>
#include <algorithm>

namespace Numerics
{
  double weightedAverage(const double & a, const double & b, const double & w_a, const double & w_b)
  {
    return (w_a * a + w_b * b) / (w_a + w_b);
  }

  double minmod(const double & a, const double & b)
  {
    const double ab = a * b;
    if (ab <= 0)
      return 0;
    else if (std::abs(a) < std::abs(b))
      return a;
    else
      return b;
  }

  double minmod(const double & a, const double & b, const double & c)
  {
    if (a > 0 && b > 0 && c > 0)
      return std::min(std::min(a, b), c);
    else if (a < 0 && b < 0 && c < 0)
      return std::max(std::max(a, b), c);
    else
      return 0;
  }

  double maxmod(const double & a, const double & b)
  {
    const double ab = a * b;
    if (ab <= 0)
      return 0;
    else if (std::abs(a) < std::abs(b))
      return b;
    else
      return a;
  }

  double computeSlopeMinmod(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R)
  {
    const double dU_L = (U_L - U) / (x_L - x);
    const double dU_R = (U_R - U) / (x_R - x);

    return minmod(dU_L, dU_R);
  }

  double computeSlopeMC(
    const double & U_L, const double & U_M, const double & U_R,
    const double & x_L, const double & x_M, const double & x_R)
  {
    const double dU_ML = (U_M - U_L) / (x_M - x_L);
    const double dU_RM = (U_R - U_M) / (x_R - x_M);
    const double dU_RL = (U_R - U_L) / (x_R - x_L);

    return minmod(dU_RL, 2 * dU_RM, 2 * dU_ML);
  }


  double computeSlopeSuperbee(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R)
  {
    const double dU_L = (U - U_L) / (x - x_L);
    const double dU_R = (U - U_R) / (x - x_R);

    const double dU1 = Numerics::minmod(dU_R, 2 * dU_L);
    const double dU2 = Numerics::minmod(2 * dU_R, dU_L);

    return maxmod(dU1, dU2);
  }

  double computeCFL(const double & dt, const double & max_wave_speed, const double & dx_min)
  {
    return max_wave_speed * dt / dx_min;
  }
}
