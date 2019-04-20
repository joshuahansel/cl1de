#ifndef DoFHandler_H
#define DoFHandler_H

#include <vector>

class DoFHandler
{
public:
  DoFHandler(const unsigned int & n_elems, const unsigned int & n_vars);

  const unsigned int & numberOfElements() const {return _n_elems;}
  const unsigned int & nVars() const {return _n_vars;}
  const unsigned int & nDoFs() const {return _n_dofs;}

  std::vector<double> getElemSolutionVector(
    const std::vector<double> & U, const unsigned int & i_elem) const;

protected:
  const unsigned int _n_elems;
  const unsigned int _n_vars;
  const unsigned int _n_dofs;
};

#endif /* DoFHandler_H */
