#ifndef ReconstructorEuler1Phase_H
#define ReconstructorEuler1Phase_H

#include <vector>

class ReconstructorEuler1Phase
{
public:
  ReconstructorEuler1Phase();

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
    std::vector<double> & rEA_R) const = 0;

protected:
  void reconstructElementSolutionGodunov(
    const std::vector<double> & rA,
    const std::vector<double> & ruA,
    const std::vector<double> & rEA,
    const std::vector<double> & A_elem,
    const std::vector<double> & A_node,
    const unsigned int & ie,
    std::vector<double> & rA_L,
    std::vector<double> & ruA_L,
    std::vector<double> & rEA_L,
    std::vector<double> & rA_R,
    std::vector<double> & ruA_R,
    std::vector<double> & rEA_R) const;
};

#endif /* ReconstructorEuler1Phase_H */
