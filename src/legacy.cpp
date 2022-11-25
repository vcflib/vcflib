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
            char* buffer = wf_aligner->cigar->operations + wf_aligner->cigar->begin_offset;
            int buf_len = wf_aligner->cigar->end_offset - wf_aligner->cigar->begin_offset;
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

/*
@@ Post-process alleles to reduce the set and normalise counts. This is the legacy version.
 */

#define ALLELE_NULL -1

void Variant::legacy_reduceAlleles(
    map<string, pair<vector<VariantAllele>, bool> > varAlleles,
    VariantCallFile &variantFile,
    Variant var,
    string parseFlag,
    bool keepInfo,
    bool keepGeno,
    bool debug)
{
    set<VariantAllele> alleles;
    // collect unique alleles
    for (auto a: varAlleles) {
        for (auto va: a.second.first) {
            if (debug) cerr << a.first << " " << va << endl;
            alleles.insert(va); // only inserts first unique allele and ignores next ones
        }
    }

    int altcount = 0;
    for (auto a: alleles) {
        if (a.ref != a.alt) {
            ++altcount;
            if (debug) cerr << altcount << "$" << a << endl;
        }
    }

    if (altcount == 1 && var.alt.size() == 1 && var.alt.front().size() == 1) { // if biallelic SNP
        cout << var << endl;
        return;
    }

    // collect variant allele indexed membership
    map<VariantAllele, vector<int> > variantAlleleIndexes; // from serialized VariantAllele to indexes
    for (auto a: varAlleles) {
        int index = var.altAlleleIndexes[a.first] + 1; // make non-relative
        for (auto va: a.second.first) {
            variantAlleleIndexes[va].push_back(index);
        }
    }

    // VariantAllele tracks pos,ref,alt. We add these counters in alleleStuff:
    struct var_info_t {
        double freq = 0;
        int count = 0;
        int in_inv = 0;
        map<string, string> info;
    };
    map<VariantAllele, var_info_t> alleleStuff;

    for (auto a: var.alt) {
        auto varalleles = varAlleles[a].first;
        bool is_inv = varAlleles[a].second;
        for (auto va: varalleles) {
            alleleStuff[va].in_inv += is_inv;
        }
    }

    bool hasAf = false;
    if (var.info.find("AF") != var.info.end()) {
        hasAf = true;
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            vector<VariantAllele>& vars = varAlleles[*a].first;
            for (vector<VariantAllele>::iterator va = vars.begin(); va != vars.end(); ++va) {
                double freq;
                try {
                    convert(var.info["AF"].at(var.altAlleleIndexes[*a]), freq);
                    alleleStuff[*va].freq += freq;
                } catch (...) {
                    cerr << "vcfallelicprimitives WARNING: AF does not have count == alts @ "
                         << var.sequenceName << ":" << var.position << endl;
                }
            }
        }
    }

    bool hasAc = false;
    if (var.info.find("AC") != var.info.end()) {
        hasAc = true;
        for (auto a: var.alt) {
            auto vars = varAlleles[a].first;
            for (auto va: vars) {
                int count;
                try {
                    convert(var.info["AC"].at(var.altAlleleIndexes[a]), count);
                    alleleStuff[va].count += count;
                } catch (...) {
                    cerr << "vcfallelicprimitives WARNING: AC does not have count == alts @ "
                         << var.sequenceName << ":" << var.position << endl;
                }
            }
        }
    }

    if (keepInfo) {
        for (map<string, vector<string> >::iterator infoit = var.info.begin();
             infoit != var.info.end(); ++infoit) {
            string key = infoit->first;
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                vector<VariantAllele>& vars = varAlleles[*a].first;
                for (vector<VariantAllele>::iterator va = vars.begin(); va != vars.end(); ++va) {
                    string val;
                    vector<string>& vals = var.info[key];
                    if (vals.size() == var.alt.size()) { // allele count for info
                        val = vals.at(var.altAlleleIndexes[*a]);
                    } else if (vals.size() == 1) { // site-wise count
                        val = vals.front();
                    } // don't handle other multiples... how would we do this without going crazy?
                    if (!val.empty()) {
                        alleleStuff[*va].info[key] = val;
                    }
                }
            }
        }
    }

    /*
      if (keepGeno) {
      for (map<string, map<string, vector<string> > >::iterator sampleit = var.samples.begin();
      sampleit != var.samples.end(); ++sampleit) {
      string& sampleName = sampleit->first;
      map<string, vector<string> >& sampleValues = var.samples[sampleName];

      }
      }
    */

    // from old allele index to a new series across the unpacked positions
    map<int, map<long unsigned int, int> > unpackedAlleleIndexes;
    map<int, bool> unpackedAlleleInversions;

    map<long unsigned int, Variant> variants;
    int varidx = 0;
    for (set<VariantAllele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        if (a->ref == a->alt) {
            // ref allele
            continue;
        }
        vector<int>& originalIndexes = variantAlleleIndexes[*a];
        string type;
        int len = 0;
        if (a->ref.size() && a->alt.size()
            && a->ref.at(0) == a->alt.at(0)) { // well-behaved indels
            if (a->ref.size() > a->alt.size()) {
                type = "del";
                len = a->ref.size() - a->alt.size();
                // special case
                // a deletion implies we should be ALLELE_NULL on this haplotype
                // until the end of the deletion
                // save the range in a new map which we'll iterate over
                for (auto i : originalIndexes) {
                    // TODO check if it should be len
                    //auto d = (*deletions)[i];
                    //d.push_back(make_pair(0, 0));
                }
            } else if (a->ref.size() < a->alt.size()) {
                len = a->alt.size() - a->ref.size();
                type = "ins";
            }
        } else {
            if (a->ref.size() == a->alt.size()) {
                len = a->ref.size();
                if (a->ref.size() == 1) {
                    type = "snp";
                } else {
                    type = "mnp";
                }
            } else {
                len = abs((int) a->ref.size() - (int) a->alt.size());
                type = "complex";
            }
        }

        if (variants.find(a->position) == variants.end()) {
            Variant newvar(variantFile);
            variants[a->position] = newvar;
        }

        Variant& v = variants[a->position]; // guaranteed to exist

        if (!parseFlag.empty()) {
            v.info[parseFlag].push_back(var.sequenceName + ":" + std::to_string(var.position));
        }
        v.quality = var.quality;
        v.filter = var.filter;
        v.infoOrderedKeys = var.infoOrderedKeys;
        if (v.id.empty()) {
            v.id = var.id + "_" + std::to_string(++varidx);
        }
        //v.format = var.format;
        vector<string> gtonlyformat;
        gtonlyformat.push_back("GT");
        v.format = gtonlyformat;
        v.info["TYPE"].push_back(type);
        v.info["LEN"].push_back(convert(len));
        v.info["INV"].push_back(convert(alleleStuff[*a].in_inv));
        if (hasAf) {
            v.info["AF"].push_back(convert(alleleStuff[*a].freq));
        }
        if (hasAc) {
            v.info["AC"].push_back(convert(alleleStuff[*a].count));
        }
        if (keepInfo) {
            for (map<string, vector<string> >::iterator infoit = var.info.begin();
                 infoit != var.info.end(); ++infoit) {
                string key = infoit->first;
                if (key != "AF" && key != "AC" && key != "TYPE" && key != "LEN") { // don't clobber previous
                    v.info[key].push_back(alleleStuff[*a].info[key]);
                }
            }
        }

        // now, keep all the other infos if we are asked to

        v.sequenceName = var.sequenceName;
        v.position = a->position; // ... by definition, this should be == if the variant was found
        if (v.ref.size() < a->ref.size()) {
            for (vector<string>::iterator va = v.alt.begin(); va != v.alt.end(); ++va) {
                *va += a->ref.substr(v.ref.size());
            }
            v.ref = a->ref;
        }
        v.alt.push_back(a->alt);

        int alleleIndex = v.alt.size();
        for (vector<int>::iterator i = originalIndexes.begin(); i != originalIndexes.end(); ++i) {
            unpackedAlleleIndexes[*i][v.position] = alleleIndex;
            //unpackedAlleleInversions[*i] = v.inv
        }
        // add null allele
        unpackedAlleleIndexes[ALLELE_NULL][v.position] = ALLELE_NULL;

    }

    // handle deletions. If ref length is larger than the WF matched
    // allele length make this a missing genotype for all individual
    // SNP/MNP calls that match the allele index and fall inside the
    // deletion.
    //
    // The idea is that when a deletion exists for a sample there is
    // no way a SNP/MNP gets called in that sample.
    for (auto a: alleles) {
        int len = 0;
        if (a.ref.size() && a.alt.size() && a.ref.at(0) == a.alt.at(0)
            && a.ref.size() > a.alt.size()) {
            len = a.ref.size() - a.alt.size();
        } else {
            continue;
        }
        assert(len > 0); // make sure we have a deletion
        // nullify all the variants inside of the deletion range by
        // walking all variants and checking the allele index
        // number. Note that this version relies on a sorted map of
        // variants[pos]. By default, a Map in C++ is sorted in
        // increasing order based on its key.
        vector<int>& originalIndexes = variantAlleleIndexes[a];
        auto begin = variants.upper_bound(a.position);
        auto end = variants.upper_bound(a.position + a.ref.size());
        for (auto i : originalIndexes) {
            for (auto x = begin; x != end; ++x) {
                unpackedAlleleIndexes[i][x->second.position] = ALLELE_NULL;
            }
        }
    }

    // genotypes
    for (vector<string>::iterator s = var.sampleNames.begin(); s != var.sampleNames.end(); ++s) {
        string& sampleName = *s;
        if (var.samples.find(sampleName) == var.samples.end()) {
            continue;
        }
        map<string, vector<string> >& sample = var.samples[sampleName];
        if (sample.find("GT") == sample.end()) {
            continue;
        }
        string& genotype = sample["GT"].front();
        vector<string> genotypeStrs = split(genotype, "|/");
        vector<int> genotypeIndexes;
        for (vector<string>::iterator s = genotypeStrs.begin(); s != genotypeStrs.end(); ++s) {
            int i;
            if (!convert(*s, i)) {
                genotypeIndexes.push_back(ALLELE_NULL);
            } else {
                genotypeIndexes.push_back(i);
            }
        }
        map<long unsigned int, vector<int> > positionIndexes;
        for (vector<int>::iterator g = genotypeIndexes.begin(); g != genotypeIndexes.end(); ++g) {
            int oldIndex = *g;
            for (map<long unsigned int, Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
                const long unsigned int& p = v->first;
                if (oldIndex == 0) { // reference
                    positionIndexes[p].push_back(0);
                } else {
                    positionIndexes[p].push_back(unpackedAlleleIndexes[oldIndex][p]);
                }
            }
        }
        for (map<long unsigned int, Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
            Variant& variant = v->second;
            vector<int>& gtints = positionIndexes[v->first];
            vector<string> gtstrs;
            for (vector<int>::iterator i = gtints.begin(); i != gtints.end(); ++i) {
                if (*i != ALLELE_NULL) {
                    gtstrs.push_back(convert(*i));
                } else {
                    gtstrs.push_back(".");
                }
            }
            string genotype = join(gtstrs, "|");
            // if we are keeping the geno info, pull it over here
            if (keepGeno) {
                variant.format = var.format;
                variant.samples[sampleName] = var.samples[sampleName];
            }
            // note that this will replace the old geno, but otherwise it is the same
            variant.samples[sampleName]["GT"].clear();
            variant.samples[sampleName]["GT"].push_back(genotype);
        }
    }
    for (auto v: variants) {
        cout << v.second << endl;
    }
}

} // namespace vcflib
