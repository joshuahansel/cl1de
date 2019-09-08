#ifndef DoFHandlerEuler2Phase_H
#define DoFHandlerEuler2Phase_H

#include "DoFHandler.h"

class DoFHandlerEuler2Phase : public DoFHandler
{
public:
  DoFHandlerEuler2Phase(const unsigned int & n_elems);

  unsigned int elemIndex_aA_liq(const unsigned int & i) const {return i * _n_vars;}
  unsigned int elemIndex_arA_liq(const unsigned int & i) const {return i * _n_vars + 1;}
  unsigned int elemIndex_aruA_liq(const unsigned int & i) const {return i * _n_vars + 2;}
  unsigned int elemIndex_arEA_liq(const unsigned int & i) const {return i * _n_vars + 3;}
  unsigned int elemIndex_arA_vap(const unsigned int & i) const {return i * _n_vars + 4;}
  unsigned int elemIndex_aruA_vap(const unsigned int & i) const {return i * _n_vars + 5;}
  unsigned int elemIndex_arEA_vap(const unsigned int & i) const {return i * _n_vars + 6;}
};

#endif /* DoFHandlerEuler2Phase_H */
