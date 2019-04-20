#include "DoFHandler.h"

DoFHandler::DoFHandler(const unsigned int & n_elems, const unsigned int & n_vars)
  : _n_elems(n_elems),
    _n_vars(n_vars),
    _n_dofs(_n_elems * _n_vars)
{
}

std::vector<double> DoFHandler::getElemSolutionVector(
  const std::vector<double> & U, const unsigned int & i_elem) const
{
  std::vector<double> U_elem(_n_vars, 0.0);
  const unsigned int base_index = i_elem * _n_vars;
  for (unsigned int i_var = 0; i_var < _n_vars; i_var++)
    U_elem[i_var] = U[base_index + i_var];

  return U_elem;
}
