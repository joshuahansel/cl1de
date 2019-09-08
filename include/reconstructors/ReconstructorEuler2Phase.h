#ifndef ReconstructorEuler2Phase_H
#define ReconstructorEuler2Phase_H

#include <vector>

class DoFHandlerEuler2Phase;

class ReconstructorEuler2Phase
{
public:
  ReconstructorEuler2Phase(const DoFHandlerEuler2Phase & dof_handler);

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

  const DoFHandlerEuler2Phase & _dof_handler;
};

#endif /* ReconstructorEuler2Phase_H */
