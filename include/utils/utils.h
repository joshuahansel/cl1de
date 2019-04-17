#ifndef utils_H
#define utils_H

#include <string>

void throwError(const std::string & message);
void throwError(const std::string & message, const std::string & prefix);
void throwInvalidStringParameterValueError(const std::string & parameter, const std::string & value);

#endif /* utils_H */
