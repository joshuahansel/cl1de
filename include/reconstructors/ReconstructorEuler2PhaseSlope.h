#ifndef ReconstructorEuler2PhaseSlope_H
#define ReconstructorEuler2PhaseSlope_H

#include "ReconstructorEuler2Phase.h"

class EOS1Phase;

class ReconstructorEuler2PhaseSlope : public ReconstructorEuler2Phase
{
public:
  ReconstructorEuler2PhaseSlope(
    const DoFHandlerEuler2Phase & dof_handler, const EOS1Phase & eos_liq, const EOS1Phase & eos_vap);

  virtual void reconstructSolution(
    const std::vector<double> & U,
    const std::vector<double> & A_elem,
    const std::vector<double> & A_node,
    const std::vector<double> & x_elem,
    const std::vector<double> & x_node,
    std::vector<double> & U_L,
    std::vector<double> & U_R) const override;

protected:
  virtual void computeSlopeVariables(
    const double & aA_liq,
    const double & arA_liq,
    const double & aruA_liq,
    const double & arEA_liq,
    const double & arA_vap,
    const double & aruA_vap,
    const double & arEA_vap,
    const double & A,
    double & w1,
    double & w2,
    double & w3,
    double & w4,
    double & w5,
    double & w6,
    double & w7) const = 0;

  virtual void computeConservativeVariables(
    const double & w1,
    const double & w2,
    const double & w3,
    const double & w4,
    const double & w5,
    const double & w6,
    const double & w7,
    const double & A,
    double & aA_liq,
    double & arA_liq,
    double & aruA_liq,
    double & arEA_liq,
    double & arA_vap,
    double & aruA_vap,
    double & arEA_vap) const = 0;

  virtual double computeSlope(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R) const = 0;

  const EOS1Phase & _eos_liq;
  const EOS1Phase & _eos_vap;
};

#endif /* ReconstructorEuler2PhaseSlope_H */
