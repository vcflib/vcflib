#pragma once

#include <vector>
#include <utility>
#include <string>

#include "var.hpp"

namespace vcflib
{
    void loadPhased(std::vector<std::pair<std::string, std::string>>& haplotypes, const genotype* pop);
}
