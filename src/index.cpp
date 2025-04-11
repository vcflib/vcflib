#include "index.hpp"

#include "split.h"

#include <vector>

#include <cstdlib>

namespace vcflib
{
void loadIndices(std::map<int, int>& index, const std::string& set)
{
	const std::vector<std::string> individuals = split(set, ",");

	for (const auto& individual : individuals)
	{
		index[std::atoi(individual.c_str())] = 1;
	}
}
}