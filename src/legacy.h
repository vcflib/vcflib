/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2024 Erik Garrison
    Copyright © 2020-2024 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#pragma once

#include "Variant.h"

#include <map>
#include <string>
#include <vector>

using namespace std;

namespace vcflib {

class VariantLegacy : public Variant {
public:

    VariantLegacy() { }

    VariantLegacy(VariantCallFile& v)
        : Variant(v)
    { }

// Legacy version of parsedAlterneates:
    map<string, vector<VariantAllele> > legacy_parsedAlternates(
	bool includePreviousBaseForIndels = false,
	bool useMNPs = false,
	bool useEntropy = false,
	float matchScore = 10.0f,
	float mismatchScore = -9.0f,
	float gapOpenPenalty = 15.0f,
	float gapExtendPenalty = 6.66f,
	float repeatGapExtendPenalty = 0.0f,
	string flankingRefLeft = "",
	string flankingRefRight = "",
	bool useWaveFront=true,
	bool debug=false);
};

} // end namespace vcflib
