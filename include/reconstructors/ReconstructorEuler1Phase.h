#ifndef ReconstructorEuler1Phase_H
#define ReconstructorEuler1Phase_H

#include <vector>

class DoFHandlerEuler1Phase;

class ReconstructorEuler1Phase
{
public:
  ReconstructorEuler1Phase(const DoFHandlerEuler1Phase & dof_handler);

  virtual void reconstructSolution(
    const std::vector<double> & U,
    const std::vector<double> & A_elem,
    const std::vector<double> & A_node,
    const std::vector<double> & x_elem,
    const std::vector<double> & x_node,
    std::vector<double> & U_L,
    std::vector<double> & U_R) const = 0;

protected:
  void reconstructElementSolutionGodunov(
    const std::vector<double> & U,
    const std::vector<double> & A_elem,
    const std::vector<double> & A_node,
    const unsigned int & ie,
    std::vector<double> & U_L,
    std::vector<double> & U_R) const;

  const DoFHandlerEuler1Phase & _dof_handler;
};

#endif /* ReconstructorEuler1Phase_H */
