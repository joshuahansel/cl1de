#include "ReconstructorEuler1Phase.h"

ReconstructorEuler1Phase::ReconstructorEuler1Phase()
{
}

void ReconstructorEuler1Phase::reconstructElementSolutionGodunov(
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
  std::vector<double> & rEA_R) const
{
  const unsigned int inL = ie;
  const unsigned int inR = ie + 1;

  const double A_ratio_L = A_node[inL] / A_elem[ie];
  rA_L[ie] = rA[ie] * A_ratio_L;
  ruA_L[ie] = ruA[ie] * A_ratio_L;
  rEA_L[ie] = rEA[ie] * A_ratio_L;

  const double A_ratio_R = A_node[inR] / A_elem[ie];
  rA_R[ie] = rA[ie] * A_ratio_R;
  ruA_R[ie] = ruA[ie] * A_ratio_R;
  rEA_R[ie] = rEA[ie] * A_ratio_R;
}
