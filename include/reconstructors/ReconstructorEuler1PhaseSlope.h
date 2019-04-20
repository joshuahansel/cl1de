#ifndef ReconstructorEuler1PhaseSlope_H
#define ReconstructorEuler1PhaseSlope_H

#include "ReconstructorEuler1Phase.h"

class EOS1Phase;

class ReconstructorEuler1PhaseSlope : public ReconstructorEuler1Phase
{
public:
  ReconstructorEuler1PhaseSlope(const DoFHandlerEuler1Phase & dof_handler, const EOS1Phase & eos);

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
    const double & rA,
    const double & ruA,
    const double & rEA,
    const double & A,
    double & w1,
    double & w2,
    double & w3) const = 0;

  virtual void computeConservativeVariables(
    const double & w1,
    const double & w2,
    const double & w3,
    const double & A,
    double & rA,
    double & ruA,
    double & rEA) const = 0;

  virtual double computeSlope(
    const double & U_L, const double & U, const double & U_R,
    const double & x_L, const double & x, const double & x_R) const = 0;

  const EOS1Phase & _eos;
};

#endif /* ReconstructorEuler1PhaseSlope_H */
