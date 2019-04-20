#include "ReconstructorEuler1Phase.h"
#include "DoFHandlerEuler1Phase.h"

ReconstructorEuler1Phase::ReconstructorEuler1Phase(const DoFHandlerEuler1Phase & dof_handler)
  : _dof_handler(dof_handler)
{
}

void ReconstructorEuler1Phase::reconstructElementSolutionGodunov(
  const std::vector<double> & U,
  const std::vector<double> & A_elem,
  const std::vector<double> & A_node,
  const unsigned int & ie,
  std::vector<double> & U_L,
  std::vector<double> & U_R) const
{
  const unsigned int inL = ie;
  const unsigned int inR = ie + 1;

  const unsigned int i_rA = _dof_handler.elemIndex_rA(ie);
  const unsigned int i_ruA = _dof_handler.elemIndex_ruA(ie);
  const unsigned int i_rEA = _dof_handler.elemIndex_rEA(ie);

  const double A_ratio_L = A_node[inL] / A_elem[ie];
  U_L[i_rA] = U[i_rA] * A_ratio_L;
  U_L[i_ruA] = U[i_ruA] * A_ratio_L;
  U_L[i_rEA] = U[i_rEA] * A_ratio_L;

  const double A_ratio_R = A_node[inR] / A_elem[ie];
  U_R[i_rA] = U[i_rA] * A_ratio_R;
  U_R[i_ruA] = U[i_ruA] * A_ratio_R;
  U_R[i_rEA] = U[i_rEA] * A_ratio_R;
}
