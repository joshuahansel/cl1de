#ifndef ReconstructorEuler1PhaseGodunov_H
#define ReconstructorEuler1PhaseGodunov_H

#include "ReconstructorEuler1Phase.h"

class ReconstructorEuler1PhaseGodunov : public ReconstructorEuler1Phase
{
public:
  ReconstructorEuler1PhaseGodunov(const DoFHandlerEuler1Phase & dof_handler);

  virtual void reconstructSolution(
    const std::vector<double> & U,
    const std::vector<double> & A_elem,
    const std::vector<double> & A_node,
    const std::vector<double> & x_elem,
    const std::vector<double> & x_node,
    std::vector<double> & U_L,
    std::vector<double> & U_R) const override;
};

#endif /* ReconstructorEuler1PhaseGodunov_H */
