#ifndef UtilsNumerics_H
#define UtilsNumerics_H

namespace Numerics
{
  double weightedAverage(const double & a, const double & b, const double & w_a, const double & w_b);

  double minmod(const double & a, const double & b);
  double minmod(const double & a, const double & b, const double & c);

  double maxmod(const double & a, const double & b);

  double computeSlopeFull(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R);
  double computeSlopeMinmod(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R);
  double computeSlopeMC(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R);
  double computeSlopeSuperbee(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R);

  double computeCFL(const double & dt, const double & max_wave_speed, const double & dx_min);
}

#endif /* UtilsNumerics_H */
