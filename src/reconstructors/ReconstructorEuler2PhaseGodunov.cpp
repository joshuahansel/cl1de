#include "ReconstructorEuler2PhaseGodunov.h"
#include "DoFHandlerEuler2Phase.h"

ReconstructorEuler2PhaseGodunov::ReconstructorEuler2PhaseGodunov(
  const DoFHandlerEuler2Phase & dof_handler)
  : ReconstructorEuler2Phase(dof_handler)
{
}

void ReconstructorEuler2PhaseGodunov::reconstructSolution(
  const std::vector<double> & U,
  const std::vector<double> & A_elem,
  const std::vector<double> & A_node,
  const std::vector<double> & /*x_elem*/,
  const std::vector<double> & /*x_node*/,
  std::vector<double> & U_L,
  std::vector<double> & U_R) const
{
  const unsigned int & n_elems = _dof_handler.numberOfElements();
  for (unsigned int ie = 0; ie < n_elems; ie++)
    reconstructElementSolutionGodunov(U, A_elem, A_node, ie, U_L, U_R);
}
