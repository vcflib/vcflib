#pragma once

#include <vector>

namespace vcflib
{
    double median(std::vector<double>& v);

    double mean(const std::vector<double>& v);
    double mean(const std::vector<int>& data);

    double variance(const std::vector <double>& data, const double mean);

    double standard_deviation(const std::vector <double>& data, const double mean);
}
