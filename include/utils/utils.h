#ifndef utils_H
#define utils_H

#include <string>

[[noreturn]] void throwError(const std::string & message);
[[noreturn]] void throwError(const std::string & message, const std::string & prefix);
[[noreturn]] void throwInvalidStringParameterValueError(const std::string & parameter, const std::string & value);

#endif /* utils_H */
