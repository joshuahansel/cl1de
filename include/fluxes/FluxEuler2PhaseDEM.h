#ifndef FluxEuler2PhaseDEM_H
#define FluxEuler2PhaseDEM_H

#include "FluxEuler2Phase.h"

class FluxEuler2PhaseDEM : public FluxEuler2Phase
{
public:
  FluxEuler2PhaseDEM(const EOS1Phase & eos_liq, const EOS1Phase & eos_vap);

  virtual void computeFlux(
    const std::vector<double> & U_L,
    const std::vector<double> & U_R,
    const double & A,
    std::vector<double> & f_L,
    std::vector<double> & f_R) const override;

protected:
  std::vector<double> computeConvectiveFlux1Phase(
    const std::vector<double> & U_L,
    const std::vector<double> & U_R,
    const EOS1Phase & eos_L,
    const EOS1Phase & eos_R,
    const double & A,
    double & u_int,
    double & p_int) const;

  void apply1PhaseFlux(
    const std::vector<double> & f_1phase,
    const double & x_times_S,
    const unsigned int & k,
    std::vector<double> & f_2phase) const;

  void apply2PhaseFlux(
    const std::vector<double> & f_2phase_add,
    const double & factor,
    std::vector<double> & f_2phase) const;

  static const unsigned int _n_phases = 2;
  const unsigned int liq = 0;
  const unsigned int vap = 1;

  static const unsigned int _n_sides = 2;
  static const unsigned int L = 0;
  static const unsigned int R = 1;
};

#endif /* FluxEuler2PhaseDEM_H */
