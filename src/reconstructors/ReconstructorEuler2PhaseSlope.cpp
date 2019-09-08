#include "ReconstructorEuler2PhaseSlope.h"
#include "DoFHandlerEuler2Phase.h"
#include "EOS1Phase.h"

ReconstructorEuler2PhaseSlope::ReconstructorEuler2PhaseSlope(
  const DoFHandlerEuler2Phase & dof_handler,
  const EOS1Phase & eos_liq,
  const EOS1Phase & eos_vap)
  : ReconstructorEuler2Phase(dof_handler),
    _eos_liq(eos_liq),
    _eos_vap(eos_vap)
{
}

void ReconstructorEuler2PhaseSlope::reconstructSolution(
  const std::vector<double> & U,
  const std::vector<double> & A_elem,
  const std::vector<double> & A_node,
  const std::vector<double> & x_elem,
  const std::vector<double> & x_node,
  std::vector<double> & U_L,
  std::vector<double> & U_R) const
{
  const unsigned int & n_elems = _dof_handler.numberOfElements();

  {
    const unsigned int ie = 0;
    reconstructElementSolutionGodunov(U, A_elem, A_node, ie, U_L, U_R);
  }
  for (unsigned int ie = 1; ie < n_elems - 1; ie++)
  {
    const unsigned int ieL = ie - 1;
    const unsigned int ieR = ie + 1;

    const unsigned int ieL_aA_liq = _dof_handler.elemIndex_aA_liq(ieL);
    const unsigned int ieL_arA_liq = _dof_handler.elemIndex_arA_liq(ieL);
    const unsigned int ieL_aruA_liq = _dof_handler.elemIndex_aruA_liq(ieL);
    const unsigned int ieL_arEA_liq = _dof_handler.elemIndex_arEA_liq(ieL);
    const unsigned int ieL_arA_vap = _dof_handler.elemIndex_arA_vap(ieL);
    const unsigned int ieL_aruA_vap = _dof_handler.elemIndex_aruA_vap(ieL);
    const unsigned int ieL_arEA_vap = _dof_handler.elemIndex_arEA_vap(ieL);

    const unsigned int ie_aA_liq = _dof_handler.elemIndex_aA_liq(ie);
    const unsigned int ie_arA_liq = _dof_handler.elemIndex_arA_liq(ie);
    const unsigned int ie_aruA_liq = _dof_handler.elemIndex_aruA_liq(ie);
    const unsigned int ie_arEA_liq = _dof_handler.elemIndex_arEA_liq(ie);
    const unsigned int ie_arA_vap = _dof_handler.elemIndex_arA_vap(ie);
    const unsigned int ie_aruA_vap = _dof_handler.elemIndex_aruA_vap(ie);
    const unsigned int ie_arEA_vap = _dof_handler.elemIndex_arEA_vap(ie);

    const unsigned int ieR_aA_liq = _dof_handler.elemIndex_aA_liq(ieR);
    const unsigned int ieR_arA_liq = _dof_handler.elemIndex_arA_liq(ieR);
    const unsigned int ieR_aruA_liq = _dof_handler.elemIndex_aruA_liq(ieR);
    const unsigned int ieR_arEA_liq = _dof_handler.elemIndex_arEA_liq(ieR);
    const unsigned int ieR_arA_vap = _dof_handler.elemIndex_arA_vap(ieR);
    const unsigned int ieR_aruA_vap = _dof_handler.elemIndex_aruA_vap(ieR);
    const unsigned int ieR_arEA_vap = _dof_handler.elemIndex_arEA_vap(ieR);

    const double & aA_liq_ieL = U[ieL_aA_liq];
    const double & arA_liq_ieL = U[ieL_arA_liq];
    const double & aruA_liq_ieL = U[ieL_aruA_liq];
    const double & arEA_liq_ieL = U[ieL_arEA_liq];
    const double & arA_vap_ieL = U[ieL_arA_vap];
    const double & aruA_vap_ieL = U[ieL_aruA_vap];
    const double & arEA_vap_ieL = U[ieL_arEA_vap];

    const double & aA_liq_ie = U[ie_aA_liq];
    const double & arA_liq_ie = U[ie_arA_liq];
    const double & aruA_liq_ie = U[ie_aruA_liq];
    const double & arEA_liq_ie = U[ie_arEA_liq];
    const double & arA_vap_ie = U[ie_arA_vap];
    const double & aruA_vap_ie = U[ie_aruA_vap];
    const double & arEA_vap_ie = U[ie_arEA_vap];

    const double & aA_liq_ieR = U[ieR_aA_liq];
    const double & arA_liq_ieR = U[ieR_arA_liq];
    const double & aruA_liq_ieR = U[ieR_aruA_liq];
    const double & arEA_liq_ieR = U[ieR_arEA_liq];
    const double & arA_vap_ieR = U[ieR_arA_vap];
    const double & aruA_vap_ieR = U[ieR_aruA_vap];
    const double & arEA_vap_ieR = U[ieR_arEA_vap];

    double w1_ieL, w2_ieL, w3_ieL, w4_ieL, w5_ieL, w6_ieL, w7_ieL;
    computeSlopeVariables(
      aA_liq_ieL, arA_liq_ieL, aruA_liq_ieL, arEA_liq_ieL, arA_vap_ieL, aruA_vap_ieL, arEA_vap_ieL, A_elem[ieL],
      w1_ieL, w2_ieL, w3_ieL, w4_ieL, w5_ieL, w6_ieL, w7_ieL);

    double w1_ie, w2_ie, w3_ie, w4_ie, w5_ie, w6_ie, w7_ie;
    computeSlopeVariables(
      aA_liq_ie, arA_liq_ie, aruA_liq_ie, arEA_liq_ie, arA_vap_ie, aruA_vap_ie, arEA_vap_ie, A_elem[ie],
      w1_ie, w2_ie, w3_ie, w4_ie, w5_ie, w6_ie, w7_ie);

    double w1_ieR, w2_ieR, w3_ieR, w4_ieR, w5_ieR, w6_ieR, w7_ieR;
    computeSlopeVariables(
      aA_liq_ieR, arA_liq_ieR, aruA_liq_ieR, arEA_liq_ieR, arA_vap_ieR, aruA_vap_ieR, arEA_vap_ieR, A_elem[ieR],
      w1_ieR, w2_ieR, w3_ieR, w4_ieR, w5_ieR, w6_ieR, w7_ieR);

    const double x_ieL = x_elem[ieL];
    const double x_ie = x_elem[ie];
    const double x_ieR = x_elem[ieR];

    const double w1_slope = computeSlope(w1_ieL, w1_ie, w1_ieR, x_ieL, x_ie, x_ieR);
    const double w2_slope = computeSlope(w2_ieL, w2_ie, w2_ieR, x_ieL, x_ie, x_ieR);
    const double w3_slope = computeSlope(w3_ieL, w3_ie, w3_ieR, x_ieL, x_ie, x_ieR);
    const double w4_slope = computeSlope(w4_ieL, w4_ie, w4_ieR, x_ieL, x_ie, x_ieR);
    const double w5_slope = computeSlope(w5_ieL, w5_ie, w5_ieR, x_ieL, x_ie, x_ieR);
    const double w6_slope = computeSlope(w6_ieL, w6_ie, w6_ieR, x_ieL, x_ie, x_ieR);
    const double w7_slope = computeSlope(w7_ieL, w7_ie, w7_ieR, x_ieL, x_ie, x_ieR);

    const unsigned int inL = ie;
    const unsigned int inR = ie + 1;

    const double dxL = x_node[inL] - x_elem[ie];
    const double dxR = x_node[inR] - x_elem[ie];

    const double w1L = w1_ie + dxL * w1_slope;
    const double w2L = w2_ie + dxL * w2_slope;
    const double w3L = w3_ie + dxL * w3_slope;
    const double w4L = w4_ie + dxL * w4_slope;
    const double w5L = w5_ie + dxL * w5_slope;
    const double w6L = w6_ie + dxL * w6_slope;
    const double w7L = w7_ie + dxL * w7_slope;

    const double w1R = w1_ie + dxR * w1_slope;
    const double w2R = w2_ie + dxR * w2_slope;
    const double w3R = w3_ie + dxR * w3_slope;
    const double w4R = w4_ie + dxR * w4_slope;
    const double w5R = w5_ie + dxR * w5_slope;
    const double w6R = w6_ie + dxR * w6_slope;
    const double w7R = w7_ie + dxR * w7_slope;

    computeConservativeVariables(
      w1L, w2L, w3L, w4L, w5L, w6L, w7L, A_node[inL],
      U_L[ie_aA_liq], U_L[ie_arA_liq], U_L[ie_aruA_liq], U_L[ie_arEA_liq], U_L[ie_arA_vap], U_L[ie_aruA_vap], U_L[ie_arEA_vap]);
    computeConservativeVariables(
      w1R, w2R, w3R, w4R, w5R, w6R, w7R, A_node[inR],
      U_R[ie_aA_liq], U_R[ie_arA_liq], U_R[ie_aruA_liq], U_R[ie_arEA_liq], U_R[ie_arA_vap], U_R[ie_aruA_vap], U_R[ie_arEA_vap]);
  }
  {
    const unsigned int ie = n_elems - 1;
    reconstructElementSolutionGodunov(U, A_elem, A_node, ie, U_L, U_R);
  }
}
