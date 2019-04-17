#include "BC.h"

BC::BC(bool is_left) : _is_left(is_left), _n(_is_left ? -1.0 : 1.0)
{
}
