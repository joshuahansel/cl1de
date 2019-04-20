#include "ReconstructorEuler1PhaseSlope.h"
#include "DoFHandlerEuler1Phase.h"
#include "EOS1Phase.h"

ReconstructorEuler1PhaseSlope::ReconstructorEuler1PhaseSlope(
  const DoFHandlerEuler1Phase & dof_handler,
  const EOS1Phase & eos)
  : ReconstructorEuler1Phase(dof_handler),
    _eos(eos)
{
}

void ReconstructorEuler1PhaseSlope::reconstructSolution(
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

    const unsigned int ieL_rA = _dof_handler.elemIndex_rA(ieL);
    const unsigned int ieL_ruA = _dof_handler.elemIndex_ruA(ieL);
    const unsigned int ieL_rEA = _dof_handler.elemIndex_rEA(ieL);

    const unsigned int ie_rA = _dof_handler.elemIndex_rA(ie);
    const unsigned int ie_ruA = _dof_handler.elemIndex_ruA(ie);
    const unsigned int ie_rEA = _dof_handler.elemIndex_rEA(ie);

    const unsigned int ieR_rA = _dof_handler.elemIndex_rA(ieR);
    const unsigned int ieR_ruA = _dof_handler.elemIndex_ruA(ieR);
    const unsigned int ieR_rEA = _dof_handler.elemIndex_rEA(ieR);

    const double & rA_ieL = U[ieL_rA];
    const double & ruA_ieL = U[ieL_ruA];
    const double & rEA_ieL = U[ieL_rEA];

    const double & rA_ie = U[ie_rA];
    const double & ruA_ie = U[ie_ruA];
    const double & rEA_ie = U[ie_rEA];

    const double & rA_ieR = U[ieR_rA];
    const double & ruA_ieR = U[ieR_ruA];
    const double & rEA_ieR = U[ieR_rEA];

    double w1_ieL, w2_ieL, w3_ieL;
    computeSlopeVariables(rA_ieL, ruA_ieL, rEA_ieL, A_elem[ieL], w1_ieL, w2_ieL, w3_ieL);

    double w1_ie, w2_ie, w3_ie;
    computeSlopeVariables(rA_ie, ruA_ie, rEA_ie, A_elem[ie], w1_ie, w2_ie, w3_ie);

    double w1_ieR, w2_ieR, w3_ieR;
    computeSlopeVariables(rA_ieR, ruA_ieR, rEA_ieR, A_elem[ieR], w1_ieR, w2_ieR, w3_ieR);

    const double x_ieL = x_elem[ieL];
    const double x_ie = x_elem[ie];
    const double x_ieR = x_elem[ieR];

    const double w1_slope = computeSlope(w1_ieL, w1_ie, w1_ieR, x_ieL, x_ie, x_ieR);
    const double w2_slope = computeSlope(w2_ieL, w2_ie, w2_ieR, x_ieL, x_ie, x_ieR);
    const double w3_slope = computeSlope(w3_ieL, w3_ie, w3_ieR, x_ieL, x_ie, x_ieR);

    const unsigned int inL = ie;
    const unsigned int inR = ie + 1;

    const double dxL = x_node[inL] - x_elem[ie];
    const double dxR = x_node[inR] - x_elem[ie];

    const double w1L = w1_ie + dxL * w1_slope;
    const double w2L = w2_ie + dxL * w2_slope;
    const double w3L = w3_ie + dxL * w3_slope;

    const double w1R = w1_ie + dxR * w1_slope;
    const double w2R = w2_ie + dxR * w2_slope;
    const double w3R = w3_ie + dxR * w3_slope;

    computeConservativeVariables(
      w1L, w2L, w3L, A_node[inL], U_L[ie_rA], U_L[ie_ruA], U_L[ie_rEA]);
    computeConservativeVariables(
      w1R, w2R, w3R, A_node[inR], U_R[ie_rA], U_R[ie_ruA], U_R[ie_rEA]);
  }
  {
    const unsigned int ie = n_elems - 1;
    reconstructElementSolutionGodunov(U, A_elem, A_node, ie, U_L, U_R);
  }
}
