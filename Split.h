// functions to split a string by a specific delimiter
#include <string>
#include <vector>
#include <sstream>
#include <string.h>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
std::vector<std::string> &split(const std::string &s, const std::string& delims, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, const std::string& delims);
