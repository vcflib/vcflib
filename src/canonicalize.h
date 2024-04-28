/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2024 Erik Garrison
    Copyright © 2020-2024 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#pragma once

#include "Variant.h"
#include <Fasta.h>

using namespace std;

namespace vcflib {

class VariantCanonical : public Variant {
public:

    VariantCanonical() { }


   /**
     * This gets set to true after canonicalize() has been called on the variant, if it succeeded.
     */
    bool canonical;

    /**
     * Get the maximum zero-based position of the reference affected by this variant.
     * Only works reliably for variants that are not SVs or for SVs that have been canonicalize()'d.
     */
    int getMaxReferencePos();

   /**
     * Convert a structural variant to the canonical VCF4.3 format using a reference.
     *   Returns true if the variant is canonicalized, false otherwise.
     *   May NOT be called twice on the same variant; it will fail an assert.
     *   Returns false for non-SVs
     *   place_seq: if true, the ref/alt fields are
     *       filled in with the corresponding sequences
     *     from the reference (and optionally insertion FASTA)
     * min_size_override: If a variant is less than this size,
     *     and it has a valid REF and ALT, consider it canonicalized
     *     even if the below conditions are not true.
     * Fully canonicalized variants (which are greater than min_size_override)
     * guarantee the following:
     *  - POS <= END and corresponds to the anchoring base for symbolic alleles
     *  - SVLEN info field is set and is positive for all variants except DELs
     *  - SVTYPE info field is set and is in {DEL, INS, INV, DUP}
     *  - END info field is set to the POS + len(REF allele) - 1 and corresponds to the final affected reference base
     *  - Insertions get an upper-case SEQ info field
     *  - REF and ALT are upper-case if filled in by this function
     *  - canonical = true;
     * TODO: CURRENTLY: canonical requires there be only one alt allele
    **/
    bool canonicalize(FastaReference& ref,
         vector<FastaReference*> insertions,
         bool place_seq = true,
         int min_size_override = 0);

    /**
     * This returns true if the variant appears able to be handled by
     * canonicalize(). It checks if it has fully specified sequence, or if it
     * has a defined SV type and length/endpoint.
     */
    bool canonicalizable();
};

} // namespace
