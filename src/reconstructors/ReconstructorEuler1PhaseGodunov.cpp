#include "ReconstructorEuler1PhaseGodunov.h"

ReconstructorEuler1PhaseGodunov::ReconstructorEuler1PhaseGodunov()
  : ReconstructorEuler1Phase()
{
}

void ReconstructorEuler1PhaseGodunov::reconstructSolution(
  const std::vector<double> & rA,
  const std::vector<double> & ruA,
  const std::vector<double> & rEA,
  const std::vector<double> & A_elem,
  const std::vector<double> & A_node,
  const std::vector<double> & /*x_elem*/,
  const std::vector<double> & /*x_node*/,
  std::vector<double> & rA_L,
  std::vector<double> & ruA_L,
  std::vector<double> & rEA_L,
  std::vector<double> & rA_R,
  std::vector<double> & ruA_R,
  std::vector<double> & rEA_R) const
{
  const unsigned int n_elems = rA.size();
  for (unsigned int ie = 0; ie < n_elems; ie++)
    reconstructElementSolutionGodunov(
      rA, ruA, rEA, A_elem, A_node, ie,
      rA_L, ruA_L, rEA_L, rA_R, ruA_R, rEA_R);
}
