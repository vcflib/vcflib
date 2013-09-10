#include "split.h"


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::string delims = std::string(1, delim);
    tokenize(s, elems, delims);
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

std::vector<std::string> &split(const std::string &s, const std::string& delims, std::vector<std::string> &elems) {
    tokenize(s, elems, delims);
    return elems;
}

std::vector<std::string> split(const std::string &s, const std::string& delims) {
    std::vector<std::string> elems;
    return split(s, delims, elems);
}
