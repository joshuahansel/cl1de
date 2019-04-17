#ifndef ReconstructorEuler1PhaseGodunov_H
#define ReconstructorEuler1PhaseGodunov_H

#include "ReconstructorEuler1Phase.h"

class ReconstructorEuler1PhaseGodunov : public ReconstructorEuler1Phase
{
public:
  ReconstructorEuler1PhaseGodunov();

  virtual void reconstructSolution(
    const std::vector<double> & rA,
    const std::vector<double> & ruA,
    const std::vector<double> & rEA,
    const std::vector<double> & A_elem,
    const std::vector<double> & A_node,
    const std::vector<double> & x_elem,
    const std::vector<double> & x_node,
    std::vector<double> & rA_L,
    std::vector<double> & ruA_L,
    std::vector<double> & rEA_L,
    std::vector<double> & rA_R,
    std::vector<double> & ruA_R,
    std::vector<double> & rEA_R) const override;
};

#endif /* ReconstructorEuler1PhaseGodunov_H */
