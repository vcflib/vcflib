/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/
#pragma once


#include <map>
#include <string>


namespace vcflib
{
    void loadIndices(std::map<int, int>& index, const std::string& set);
}
