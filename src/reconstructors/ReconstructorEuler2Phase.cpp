#include "ReconstructorEuler2Phase.h"
#include "DoFHandlerEuler2Phase.h"

ReconstructorEuler2Phase::ReconstructorEuler2Phase(const DoFHandlerEuler2Phase & dof_handler)
  : _dof_handler(dof_handler)
{
}

void ReconstructorEuler2Phase::reconstructElementSolutionGodunov(
  const std::vector<double> & U,
  const std::vector<double> & A_elem,
  const std::vector<double> & A_node,
  const unsigned int & ie,
  std::vector<double> & U_L,
  std::vector<double> & U_R) const
{
  const unsigned int inL = ie;
  const unsigned int inR = ie + 1;

  const unsigned int i_aA_liq = _dof_handler.elemIndex_aA_liq(ie);
  const unsigned int i_arA_liq = _dof_handler.elemIndex_arA_liq(ie);
  const unsigned int i_aruA_liq = _dof_handler.elemIndex_aruA_liq(ie);
  const unsigned int i_arEA_liq = _dof_handler.elemIndex_arEA_liq(ie);
  const unsigned int i_arA_vap = _dof_handler.elemIndex_arA_vap(ie);
  const unsigned int i_aruA_vap = _dof_handler.elemIndex_aruA_vap(ie);
  const unsigned int i_arEA_vap = _dof_handler.elemIndex_arEA_vap(ie);

  const double A_ratio_L = A_node[inL] / A_elem[ie];
  U_L[i_aA_liq] = U[i_aA_liq] * A_ratio_L;
  U_L[i_arA_liq] = U[i_arA_liq] * A_ratio_L;
  U_L[i_aruA_liq] = U[i_aruA_liq] * A_ratio_L;
  U_L[i_arEA_liq] = U[i_arEA_liq] * A_ratio_L;
  U_L[i_arA_vap] = U[i_arA_vap] * A_ratio_L;
  U_L[i_aruA_vap] = U[i_aruA_vap] * A_ratio_L;
  U_L[i_arEA_vap] = U[i_arEA_vap] * A_ratio_L;

  const double A_ratio_R = A_node[inR] / A_elem[ie];
  U_R[i_aA_liq] = U[i_aA_liq] * A_ratio_R;
  U_R[i_arA_liq] = U[i_arA_liq] * A_ratio_R;
  U_R[i_aruA_liq] = U[i_aruA_liq] * A_ratio_R;
  U_R[i_arEA_liq] = U[i_arEA_liq] * A_ratio_R;
  U_R[i_arA_vap] = U[i_arA_vap] * A_ratio_R;
  U_R[i_aruA_vap] = U[i_aruA_vap] * A_ratio_R;
  U_R[i_arEA_vap] = U[i_arEA_vap] * A_ratio_R;
}
