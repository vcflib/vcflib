/*
    vcflib C++ library for parsing and manipulating VCF files. This file contains
    legacy material that will be phased out.

    Copyright © 2010-2022 Erik Garrison
    Copyright © 2020-2022 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include <utility>
#include "Variant.h"
#include "allele.hpp"
#include "cigar.hpp"

namespace vcflib {

// parsedAlternates returns a hash of 'ref' and a vector of alts. A
// single record may be split into multiple records with new
// 'refs'. In this function Smith-Waterman is used with padding on
// both sides of a ref and each alt. The SW method quadratic in nature
// and painful with long sequences. Recently WFA is introduced with
// runs in linear time.
//
// Returns map of [REF,ALTs] with attached VariantAllele records

map<string, vector<VariantAllele> > Variant::legacy_parsedAlternates(
    bool includePreviousBaseForIndels,
    bool useMNPs,
    bool useEntropy,
    float matchScore,

    float mismatchScore,
    float gapOpenPenalty,
    float gapExtendPenalty,
    float repeatGapExtendPenalty,
    string flankingRefLeft,
    string flankingRefRight,
    bool useWaveFront,
    bool debug) {

    bool OLDPADDING = true; // we want to use the first bp for alignment

    map<string, vector<VariantAllele> > variantAlleles; // return type

    if (isSymbolicSV()){
        // Don't ever align symbolic SVs. It just wrecks things.
        return this->flatAlternates();
    }
    // add the reference allele
    variantAlleles[ref].push_back(VariantAllele(ref, ref, position));

    // single SNP case, no ambiguity possible, no need to spend a lot of
    // compute aligning ref and alt fields
    if (alt.size() == 1 && ref.size() == 1 && alt.front().size() == 1) {
        variantAlleles[alt.front()].push_back(VariantAllele(ref, alt.front(), position));
        return variantAlleles;
    }

    // Padding is used to ensure a stable alignment of the alternates to the reference
    // without having to go back and look at the full reference sequence.
    // Dynamically determine optimum padding length; if the ref is
    // larger than 10 we take a larger size for SW (the largest of ref
    // and alt, see below)
    string lpadding, rpadding;
    int paddingLen = -1;
    char anchorChar = 'Q'; // overwrites first pos of sequence

    if (flankingRefLeft.empty() && flankingRefRight.empty()) {
        paddingLen = (useWaveFront ? 10 : max(10, (int) (ref.size())));
        for (auto a: alt) {
            paddingLen = max(paddingLen, (int) (a.size()));
        }
        char padChar = 'Z';
        rpadding = string(paddingLen, padChar); // repeat padChar
        lpadding = rpadding;
    }
    else {
        rpadding = flankingRefRight;
        lpadding = flankingRefLeft;
        paddingLen = lpadding.size();
    }

    // this 'anchored' string is done for stability
    // the assumption is that there should be a positional match in the first base
    // this is true for VCF 4.1, and standard best practices
    // using the anchor char ensures this without other kinds of realignment
    string reference_M = lpadding + anchorChar + ref + rpadding;
    if (OLDPADDING) {
        reference_M = lpadding + ref + rpadding;
        reference_M[paddingLen] = anchorChar; // patch sequence with anchor
    }

    if (debug) cerr << "====>" << reference_M << endl;
    // ACCCCCACCCCCACC
    // padded ref ZZZZZZZZZZZZZZZQACCCCCACCCCCACCZZZZZZZZZZZZZZZ

    for (auto alternate: alt) { // iterate ALT strings
        unsigned int referencePos;
        string alternateQuery_M;
        if (OLDPADDING) {
          alternateQuery_M = lpadding + alternate + rpadding;
          alternateQuery_M[paddingLen] = anchorChar; // patch sequence with anchor
        }
        else
            alternateQuery_M = lpadding + anchorChar + alternate + rpadding;

        if (debug) cerr << alternate << " => " << alternateQuery_M << endl;

        /* ['ACC', 'AC', 'ACCCCCACCCCCAC', 'ACCCCCACC', 'ACA']
           1: ACC => ZZZZZZZZZZZZZZZQCCZZZZZZZZZZZZZZZ
           1: AC => ZZZZZZZZZZZZZZZQCZZZZZZZZZZZZZZZ
           1: ACCCCCACCCCCAC => ZZZZZZZZZZZZZZZQCCCCCACCCCCACZZZZZZZZZZZZZZZ
           1: ACCCCCACC => ZZZZZZZZZZZZZZZQCCCCCACCZZZZZZZZZZZZZZZ
           1: ACA => ZZZZZZZZZZZZZZZQCAZZZZZZZZZZZZZZZ
        */

        string cigar;
        vector<pair<int, char> > cigarData;

        if (useWaveFront)
        {
            /*
             * WFA2-lib
             */
            /*
            // the C++ WFA2-lib interface is not yet stable due to heuristic initialization issues
            WFAlignerGapAffine2Pieces aligner(19,39,3,81,1,WFAligner::Alignment,WFAligner::MemoryHigh);
            aligner.alignEnd2End(reference_M.c_str(), reference_M.size(), alternateQuery_M.c_str(), alternateQuery_M.size());
            cigar = aligner.getAlignmentCigar();
            */
            auto attributes = wavefront_aligner_attr_default;
            attributes.memory_mode = wavefront_memory_low;
            attributes.distance_metric = gap_affine_2p;
            attributes.affine2p_penalties.match = 0;
            attributes.affine2p_penalties.mismatch = 19;
            attributes.affine2p_penalties.gap_opening1 = 39;
            attributes.affine2p_penalties.gap_extension1 = 3;
            attributes.affine2p_penalties.gap_opening2 = 81;
            attributes.affine2p_penalties.gap_extension2 = 1;
            attributes.alignment_scope = compute_alignment;
            auto wf_aligner = wavefront_aligner_new(&attributes);
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
            char* buffer = wf_aligner->cigar.operations + wf_aligner->cigar.begin_offset;
            int buf_len = wf_aligner->cigar.end_offset - wf_aligner->cigar.begin_offset;
            // Create string and return
            cigar = std::string(buffer,buf_len);
            if (debug) cerr << "Have a CIGAR " << cigar << endl;
            wavefront_aligner_delete(wf_aligner);
            //if (debug)
            //    cerr << "WFA output [" << cigar << "]" << endl;
            if (cigar == "") {
                if (debug) {
                    cerr << "Warning: skipping input with WF because there is no CIGAR!" << endl;
                }
                cerr << ">fail.pattern" << endl
                     << reference_M << endl
                     << ">fail.query" << endl
                     << alternateQuery_M << endl;
                variantAlleles[alt.front()].push_back(VariantAllele(ref, alt.front(), position));
                return variantAlleles;
            }
            cigarData = splitUnpackedCigar(cigar);
        }
        else {
            CSmithWatermanGotoh sw(matchScore,
                                   mismatchScore,
                                   gapOpenPenalty,
                                   gapExtendPenalty);
            if (useEntropy) sw.EnableEntropyGapPenalty(1);
            if (repeatGapExtendPenalty != 0){
                sw.EnableRepeatGapExtensionPenalty(repeatGapExtendPenalty);
            }
            sw.Align(referencePos, cigar, reference_M, alternateQuery_M);
            if (debug) cerr << "SW CIGAR " << cigar << endl;
            cigarData = splitCigar(cigar);
        }

        if (debug)
          cerr << (useWaveFront ? "WF=" : "SW=") << referencePos << ":" << cigar << ":" << reference_M << "," << alternateQuery_M << endl;

        if (cigarData.size() == 0) {
            cerr << "Algorithm error: CIGAR <" << cigar << "> is empty for "
                 << "ref " << reference_M << ","
                 << "allele " << alternateQuery_M << endl;

            exit(1);
        }

        // Check for matched padding (ZZZs)
        if (cigarData.front().second != 'M'
            || cigarData.back().second != 'M'
            || cigarData.front().first < paddingLen
            || cigarData.back().first < paddingLen) {
            cerr << "parsedAlternates: alignment does not start or end with match over padded  sequence" << endl;
            cerr << "pos: " << position << endl;
            cerr << "cigar: " << cigar << endl;
            cerr << "cigardata: " << joinCigar(cigarData) << endl;
            cerr << "ref:   " << reference_M << endl;
            cerr << "allele:" << alternateQuery_M << endl;
            exit(1);
        } else {
            // Remove the padding and anchor character
            cigarData.front().first -= paddingLen;
            cigarData.back().first -= paddingLen;
            if (!OLDPADDING) cigarData.front().first -= 1;
        }
        cigar = joinCigar(cigarData);

        if (debug)
          cerr << referencePos << ":" << cigar << ":" << reference_M << "," << alternateQuery_M << endl;

        // Walk the CIGAR for one alternate and build up variantAlleles
        vector<VariantAllele> &variants = variantAlleles[alternate];
        if (!variants.empty()) {
            cerr << "Found duplicate allele and ignoring it! " << alternate << endl;
            cerr << "pos: " << position << endl;
            continue;
        }
        int altpos = 0;
        int refpos = 0;

        for (auto e: cigarData) {
            int  mlen  = e.first;  // CIGAR matchlen
            char mtype = e.second; // CIGAR matchtype

            switch (mtype) {
            case 'I': // CIGAR INSERT
                if (includePreviousBaseForIndels) {
                    if (!variants.empty() &&
                        variants.back().ref != variants.back().alt) {
                        VariantAllele a =
                            VariantAllele("",
                                          alternate.substr(altpos, mlen),
                                          refpos + position);
                        variants.back() = variants.back() + a;
                    } else {
                        VariantAllele a =
                            VariantAllele(ref.substr(refpos - 1, 1),
                                          alternate.substr(altpos - 1, mlen + 1),
                                          refpos + position - 1);
                        variants.push_back(a);
                    }
                } else {
                    // inject insertion from allele
                    variants.push_back(VariantAllele("",
                                                     alternate.substr(altpos, mlen),
                                                     refpos + position));
                }
                altpos += mlen;
                break;
            case 'D': // CIGAR DELETE
                if (includePreviousBaseForIndels) {
                    if (!variants.empty() &&
                        variants.back().ref != variants.back().alt) {
                        VariantAllele a
                            = VariantAllele(ref.substr(refpos, mlen)
                                            , "", refpos + position);
                        variants.back() = variants.back() + a;
                    } else {
                        VariantAllele a
                            = VariantAllele(ref.substr(refpos - 1, mlen + 1),
                                            alternate.substr(altpos - 1, 1),
                                            refpos + position - 1);
                        variants.push_back(a);
                    }
                } else {
                    // inject deletion from ref
                    variants.push_back(VariantAllele(ref.substr(refpos, mlen),
                                                     "", refpos + position));
                }
                refpos += mlen;
                break;
            case 'M': // CIGAR match and variant
            case 'X':
            {
                for (int i = 0; i < mlen; ++i) {
                    VariantAllele a
                        = VariantAllele(ref.substr(refpos + i, 1),
                                        alternate.substr(altpos + i, 1),
                                        (refpos + i + position));
                    if (useMNPs) {
                        bool is_first = variants.empty();
                        auto &tail = variants.back(); // note that with is_first tail is undefined
                        auto &ref = tail.ref;
                        auto &alt = tail.alt;

                        if (!is_first && ref.size() == alt.size() && ref != alt) {
                            // elongate last variant as MNP
                            tail = tail + a;
                        } else {
                            // Inject variants as a single nucleotides
                            variants.push_back(a);
                        }
                    }
                }
            }
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
    }
    return variantAlleles;
}

} // namespace vcflib
