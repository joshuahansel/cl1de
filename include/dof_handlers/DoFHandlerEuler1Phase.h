#ifndef DoFHandlerEuler1Phase_H
#define DoFHandlerEuler1Phase_H

#include "DoFHandler.h"

class DoFHandlerEuler1Phase : public DoFHandler
{
public:
  DoFHandlerEuler1Phase(const unsigned int & n_elems);

  unsigned int elemIndex_rA(const unsigned int & i) const {return i * _n_vars;}
  unsigned int elemIndex_ruA(const unsigned int & i) const {return i * _n_vars + 1;}
  unsigned int elemIndex_rEA(const unsigned int & i) const {return i * _n_vars + 2;}
};

#endif /* DoFHandlerEuler1Phase_H */
