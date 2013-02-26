#ifndef __CONVERT_H
#define __CONVERT_H

#include <sstream>

// converts the string into the specified type, setting r to the converted
// value and returning true/false on success or failure
template<typename T>
bool inline convert(const std::string& s, T& r) {
    std::istringstream iss(s);
    iss >> r;
    return iss.eof() ? true : false;
}

template<typename T>
std::string inline convert(const T& r) {
    std::ostringstream oss;
    oss << r;
    return oss.str();
}

#endif
