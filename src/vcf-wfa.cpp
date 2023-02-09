/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2023 Erik Garrison
    Copyright © 2020-2023 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include "vcf-wfa.h"

namespace vcflib {

// parsedAlternates returns a hash of 'ref' and a vector of alts. A
// single record may be split into multiple records with new
// 'refs'. In this function Smith-Waterman is used with padding on
// both sides of a ref and each alt. The SW method quadratic in nature
// and painful with long sequences. Recently WFA is introduced with
// runs in linear time.
//
// Returns map of [REF,ALTs] with attached VariantAllele records


map<string, pair<vector<VariantAllele>,bool> > WfaVariant::parsedAlternates(
    bool includePreviousBaseForIndels,
    bool useMNPs,
    bool useEntropy,
    string flankingRefLeft,
    string flankingRefRight,
    wavefront_aligner_attr_t* wfaParams,
    int invKmerLen,
    int invMinLen,
    int threads,
    bool debug) {

    // return type is a hash of allele containing a list of variants and a
    // boolean value signifying an inversion of the variant
    map<string, pair<vector<VariantAllele>,bool> > variantAlleles;

    if (isSymbolicSV()){
        // Don't ever align symbolic SVs. It just wrecks things.
        for (auto& f : this->flatAlternates()) {
            variantAlleles[f.first] = make_pair(f.second, false);
        }
        return variantAlleles;
    }
    // add the reference allele
    variantAlleles[ref].first.push_back(VariantAllele(ref, ref, position));

    // single SNP case, no ambiguity possible, no need to spend a lot of
    // compute aligning ref and alt fields
    if (alt.size() == 1 && ref.size() == 1 && alt.front().size() == 1) {
        variantAlleles[alt.front()].first.push_back(VariantAllele(ref, alt.front(), position));
        return variantAlleles;
    }

    const string& reference_M = ref;

    std::vector<rkmh::hash_t> ref_sketch_fwd;
    std::vector<rkmh::hash_t> ref_sketch_rev;
    if (invKmerLen && ref.size() >= invKmerLen
        && ref.size() >= invMinLen) {
        ref_sketch_fwd = rkmh::hash_sequence(
            ref.c_str(), ref.size(), invKmerLen, ref.size()-invKmerLen+1);
        string fer = reverse_complement(ref);
        ref_sketch_rev = rkmh::hash_sequence(
            fer.c_str(), fer.size(), invKmerLen, fer.size()-invKmerLen+1);

    }

    /*
    for (auto a: alt) { // iterate ALT strings
        // unsigned int referencePos;
        string& alternate = a;
        variantAlleles[alternate]; // create the slot
    }
    */

// #pragma omp parallel for
    for (auto a: alt) { // iterate ALT strings
        //for (uint64_t idx = 0; idx < alt.size(); ++idx) {
        //auto& a = alt[idx];
        // unsigned int referencePos;
        string alternate = a;
        pair<vector<VariantAllele>, bool>& _v = variantAlleles[alternate];
        bool is_inv = false;
        // get the alt/ref mapping in inversion space
        if (invKmerLen && alternate.size() >= invKmerLen
            && alternate.size() >= invMinLen
            && ref.size() >= invKmerLen
            && ref.size() >= invMinLen
            && max((int)ref.size(), (int)alt.size()) > invMinLen) {
            // check if it's more likely for us to align as an inversion
            auto alt_sketch = rkmh::hash_sequence(
                a.c_str(), a.size(), invKmerLen, a.size()-invKmerLen+1);
            if (rkmh::compare(alt_sketch, ref_sketch_fwd, invKmerLen)
                > rkmh::compare(alt_sketch, ref_sketch_rev, invKmerLen)) {
                is_inv = true;
                // flip the alt
                string alternate_rev = reverse_complement(alternate);
                // cout << "alt: " << alternate_rev << endl;
                alternate = alternate_rev;
            }
            variantAlleles[alternate]; // create slot
            variantAlleles[alternate].second = is_inv;
            /*
            cerr << "comparing "
                 << " vs ref fwd " << rkmh::compare(alt_sketch, ref_sketch_fwd, invKmerLen)
                 << " vs rev " << rkmh::compare(alt_sketch, ref_sketch_rev, invKmerLen) << endl;
            */
        }

        const string& alternateQuery_M = alternate;

        string cigar;
        vector<pair<int, char> > cigarData;

        /*
         * WFA2-lib
         */
        /*
        // the C++ WFA2-lib interface is not yet stable due to heuristic mode initialization issues
        WFAlignerGapAffine2Pieces aligner(19,39,3,81,1,WFAligner::Alignment,WFAligner::MemoryHigh);
        aligner.alignEnd2End(reference_M.c_str(), reference_M.size(), alternateQuery_M.c_str(), alternateQuery_M.size());
        cigar = aligner.getAlignmentCigar();
        */
        auto wfp = *wfaParams;
        size_t max_len = max(reference_M.size(), alternateQuery_M.size());
        if (max_len > 1000) {
            wfp.memory_mode = wavefront_memory_ultralow;
        } else {
            wfp.memory_mode = wavefront_memory_high;
        }
        auto wf_aligner = wavefront_aligner_new(&wfp);

        wavefront_aligner_set_max_num_threads(wf_aligner,threads);
        wavefront_aligner_set_heuristic_none(wf_aligner);
        wavefront_aligner_set_alignment_end_to_end(wf_aligner);
        wavefront_align(wf_aligner,
                        reference_M.c_str(), reference_M.size(),
                        alternateQuery_M.c_str(), alternateQuery_M.size());
        //cerr << "WFA_input\t" << ">
        /*
          cigar_print_pretty(stderr,
          reference_M.c_str(), reference_M.size(),
          alternateQuery_M.c_str(), alternateQuery_M.size(),
          &wf_aligner->cigar,wf_aligner->mm_allocator);
        */
        // Fetch CIGAR
        char* buffer = wf_aligner->cigar->operations + wf_aligner->cigar->begin_offset;
        int buf_len = wf_aligner->cigar->end_offset - wf_aligner->cigar->begin_offset;
        // Create string and return
        cigar = std::string(buffer,buf_len);
        wavefront_aligner_delete(wf_aligner);
        if (cigar == "") {
            if (debug) {
                cerr << "ERROR: stopped WF because there is no CIGAR!" << endl;
            }
            cerr << ">fail.pattern" << endl
                 << reference_M << endl
                 << ">fail.query" << endl
                 << alternateQuery_M << endl;
            variantAlleles[alt.front()].first.push_back(VariantAllele(ref, alt.front(), position));
            exit(1);
        }
        cigarData = splitUnpackedCigar(cigar);

        if (debug)
            cerr << "biWF= " << position << ":" << cigar << "/" << joinCigar(cigarData) << ":" << reference_M << "," << alternateQuery_M << endl;

        if (cigarData.size() == 0) {
            cerr << "Algorithm error: CIGAR <" << cigar << "> is empty for "
                 << "ref " << reference_M << ","
                 << "allele " << alternateQuery_M << endl;

            exit(1);
        }

        // now left align!
        //
        // TODO: currently broken! it seems to mess up our indel alleles (they become length 0 on ref or alt)
        //stablyLeftAlign(alternateQuery_M, reference_M, cigarData, 5, true);

        if (debug) {
            cigar = joinCigar(cigarData);
            cerr << position << ":" << cigar << ":" << reference_M << "," << alternateQuery_M << endl;
        }

        // Walk the CIGAR for one alternate and build up variantAlleles
        vector<VariantAllele> &variants = variantAlleles[alternate].first;
        int altpos = 0;
        int refpos = 0;

        for (auto e: cigarData) {
            int  mlen  = e.first;  // CIGAR matchlen
            char mtype = e.second; // CIGAR matchtype

            switch (mtype) {
            case 'I': // CIGAR INSERT
            {
                auto allele = VariantAllele("",
                                            alternate.substr(altpos, mlen),
                                            refpos + position);
                // inject insertion from allele
                if (variants.size() && variants.back().is_pure_indel()) {
                    variants.back() = variants.back() + allele;
                } else {
                    variants.push_back(allele);
                }
                altpos += mlen;
            }
            break;
            case 'D': // CIGAR DELETE
            {
                auto allele = VariantAllele(ref.substr(refpos, mlen),
                                            "", refpos + position);
                // inject deletion from ref
                if (variants.size() && variants.back().is_pure_indel()) {
                    variants.back() = variants.back() + allele;
                } else {
                    variants.push_back(allele);
                }
                refpos += mlen;
            }
            break;
            case 'M': // CIGAR match and variant
            case 'X':
                variants.push_back(VariantAllele(ref.substr(refpos, mlen),
                                                 alternate.substr(altpos, mlen),
                                                 refpos + position));
                refpos += mlen;
                altpos += mlen;
                break;
            case 'S': // S operations specify segments at the start
                      // and/or end of the query that do not appear in
                      // a local alignment
                refpos += mlen;
                altpos += mlen;
                break;

            default: // Ignore N, H
                cerr << "Hit unexpected CIGAR " << mtype << " in " << endl;
                cerr << "pos: " << position << endl;
                cerr << "cigar: " << cigar << endl;
                cerr << "ref:   " << reference_M << endl;
                cerr << "allele:" << alternateQuery_M << endl;
                if (debug) exit(1);
                break;
            } // switch mtype
        }
        if (includePreviousBaseForIndels) {
            for (uint64_t i = 0; i < variants.size(); ++i) {
                auto& v = variants[i];
                if (v.is_pure_indel()) {
                    if (i == 0) { // special case, we're at the beginning
                        // do we have a next allele?
                        if (i == variants.size()-1) {
                            // if not-- panic
                            cerr << "allele base fail: no subsequent alleles to indel" << endl;
                            cerr << v << endl;
                            cerr << "variant at " << position << endl;
                            exit(1);
                        }
                        // else
                        // while the next allele is an indel
                        auto j = i+1;
                        while (j < variants.size()
                               // stop when we have sequence in both ref or alt
                               && v.is_pure_indel()) {
                            auto q = variants[j++];
                            // merge with us
                            shift_mid_left(v, q);
                        }
                        // panic if we reach the end of the variants
                        if (j == variants.size() && v.is_pure_indel()) {
                            cerr << "allele base fail: can't get an additional base next" << endl;
                            cerr << v << endl;
                            cerr << "variant at " << position << endl;
                            exit(1);
                        }
                    } else {
                        // while the next allele is an indel
                        int j = i-1; // i is guaranteed >=1
                        while (j >= 0
                               // stop when we have sequence in both ref or alt
                               && v.is_pure_indel()) {
                            auto q = variants[j--];
                            // move the middle base between the alleles
                            // to be at the start of v
                            // or merge if both are pure indels
                            shift_mid_right(q, v);
                        }
                        // panic if we reach the end of the variants
                        if (j <= 0 && v.is_pure_indel()) {
                            cerr << "allele base fail: can't get an additional base prev" << endl;
                            cerr << v << endl;
                            cerr << "variant at " << position << endl;
                            exit(1);
                        }
                        // while the previous allele is an indel
                        // merge with us
                        // stop when we have sequence in both ref or alt
                        // if allele we reach is a positional match
                        // panic if we reach the start of the variants

                    }
                    // try to get the previous reference or alternate base
                    // continue until the
                    // we first check if we have an indel just before us
                    // if so, we need to merge this into the same allele and go further back
                }
            }
        }
    }

    return variantAlleles;
}

} // namespace
