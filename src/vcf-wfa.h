/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2023 Erik Garrison
    Copyright © 2020-2023 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#pragma once

#include "wavefront/wfa.hpp"

namespace vcflib {

class WfaVariant : public Variant
{

public:
    WfaVariant(VariantCallFile& v) : Variant(v)
      { };

    map<string, pair<vector<VariantAllele>,bool> > parsedAlternates(
           bool includePreviousBaseForIndels = false,
           bool useMNPs = false,
           bool useEntropy = false,
           string flankingRefLeft = "",
           string flankingRefRight = "",
           wavefront_aligner_attr_t* wfaParams = NULL,
           int invKmerLen = 17,
           int invMinLen = 1000,
           int threads = 1,
           bool debug = false);
};

}
