#include <iostream>

#include "utils.h"

void throwError(const std::string & message)
{
  std::cerr << "\033[31mERROR:\n  "<< message << "\033[0m" << std::endl << std::flush;

  std::abort();
}

void throwError(const std::string & message, const std::string & prefix)
{
  throwError(prefix + ":\n    " + message);
}

void throwInvalidStringParameterValueError(const std::string & parameter, const std::string & value)
{
  throwError("Invalid value '" + value + "' for parameter '" + parameter + "'.");
}
