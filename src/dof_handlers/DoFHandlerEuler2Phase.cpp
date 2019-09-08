#include "DoFHandlerEuler2Phase.h"

DoFHandlerEuler2Phase::DoFHandlerEuler2Phase(const unsigned int & n_elems)
  : DoFHandler(n_elems, 7)
{
}
