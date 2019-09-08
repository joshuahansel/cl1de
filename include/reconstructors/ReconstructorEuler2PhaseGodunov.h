#ifndef ReconstructorEuler2PhaseGodunov_H
#define ReconstructorEuler2PhaseGodunov_H

#include "ReconstructorEuler2Phase.h"

class ReconstructorEuler2PhaseGodunov : public ReconstructorEuler2Phase
{
public:
  ReconstructorEuler2PhaseGodunov(const DoFHandlerEuler2Phase & dof_handler);

  virtual void reconstructSolution(
    const std::vector<double> & U,
    const std::vector<double> & A_elem,
    const std::vector<double> & A_node,
    const std::vector<double> & x_elem,
    const std::vector<double> & x_node,
    std::vector<double> & U_L,
    std::vector<double> & U_R) const override;
};

#endif /* ReconstructorEuler2PhaseGodunov_H */
