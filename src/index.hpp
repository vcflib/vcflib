/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/
#pragma once

#include "split.h"

#include <map>
#include <string>
#include <vector>

#include <cstdlib>

namespace vcflib
{
	inline void loadIndices(std::map<int, int>& index, const std::string& set)
	{
		const std::vector<std::string> individuals = split(set, ",");

		for (const auto& individual : individuals)
		{
			index[std::atoi(individual.c_str())] = 1;
		}
	}
}
