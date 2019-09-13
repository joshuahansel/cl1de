#ifndef FluxEuler2PhaseHLLC_H
#define FluxEuler2PhaseHLLC_H

#include "FluxEuler2Phase.h"

class FluxEuler2PhaseHLLC : public FluxEuler2Phase
{
public:
  FluxEuler2PhaseHLLC(const EOS1Phase & eos_liq, const EOS1Phase & eos_vap);

  virtual void computeFlux(
    const std::vector<double> & U_L,
    const std::vector<double> & U_R,
    const double & A,
    std::vector<double> & f_L,
    std::vector<double> & f_R) const override;

  std::vector<std::vector<double>> computeConservativeFlux(
    const std::vector<double> & U_L,
    const std::vector<double> & U_R,
    const double & A) const;

  std::vector<unsigned int> getLastRegionIndices() const {return _last_region_indices;}

protected:
  std::vector<double> solutionSupersonic(const std::vector<double> & W, const double & A) const;

  std::vector<double> solutionSubsonic(
    const std::vector<double> & W,
    const double & S,
    const double & SM,
    const double & A) const;

  std::vector<double> solutionSubsonicInterfacialLeft(
    const std::vector<double> & WL,
    const std::vector<double> & WR,
    const double & SL,
    const double & SR,
    const double & SM,
    const double & p_int,
    const double & A) const;

  void applyRankineHugoniot(
    std::vector<double> & flux,
    const double & S,
    const std::vector<double> & U1,
    const std::vector<double> & U2) const;

  mutable std::vector<unsigned int> _last_region_indices;
  mutable std::vector<std::vector<double>> _a;
  mutable double _u_int;
  mutable double _p_int;

  static const unsigned int _n_phases = 2;
  static const unsigned int liq = 0;
  static const unsigned int vap = 1;

  static const unsigned int _n_sides = 2;
  static const unsigned int L = 0;
  static const unsigned int R = 1;
};

#endif /* FluxEuler2PhaseHLLC_H */
