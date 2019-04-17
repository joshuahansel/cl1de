#include "ReconstructorEuler1PhaseSlope.h"
#include "EOS1Phase.h"

ReconstructorEuler1PhaseSlope::ReconstructorEuler1PhaseSlope(const EOS1Phase & eos)
  : ReconstructorEuler1Phase(),
    _eos(eos)
{
}

void ReconstructorEuler1PhaseSlope::reconstructSolution(
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
  std::vector<double> & rEA_R) const
{
  const unsigned int n_elems = rA.size();

  {
    const unsigned int ie = 0;
    reconstructElementSolutionGodunov(
      rA, ruA, rEA, A_elem, A_node, ie,
      rA_L, ruA_L, rEA_L, rA_R, ruA_R, rEA_R);
  }
  for (unsigned int ie = 1; ie < n_elems - 1; ie++)
  {
    const unsigned int ieL = ie - 1;
    const unsigned int ieR = ie + 1;

    double w1_ieL, w2_ieL, w3_ieL;
    computeSlopeVariables(rA[ieL], ruA[ieL], rEA[ieL], A_elem[ieL], w1_ieL, w2_ieL, w3_ieL);

    double w1_ie, w2_ie, w3_ie;
    computeSlopeVariables(rA[ie], ruA[ie], rEA[ie], A_elem[ie], w1_ie, w2_ie, w3_ie);

    double w1_ieR, w2_ieR, w3_ieR;
    computeSlopeVariables(rA[ieR], ruA[ieR], rEA[ieR], A_elem[ieR], w1_ieR, w2_ieR, w3_ieR);

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
      w1L, w2L, w3L, A_node[inL], rA_L[ie], ruA_L[ie], rEA_L[ie]);
    computeConservativeVariables(
      w1R, w2R, w3R, A_node[inR], rA_R[ie], ruA_R[ie], rEA_R[ie]);
  }
  {
    const unsigned int ie = n_elems - 1;
    reconstructElementSolutionGodunov(
      rA, ruA, rEA, A_elem, A_node, ie,
      rA_L, ruA_L, rEA_L, rA_R, ruA_R, rEA_R);
  }
}
