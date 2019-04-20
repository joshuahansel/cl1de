#include "DoFHandlerEuler1Phase.h"

DoFHandlerEuler1Phase::DoFHandlerEuler1Phase(const unsigned int & n_elems)
  : DoFHandler(n_elems, 3)
{
}
