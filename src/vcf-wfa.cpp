/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2024= Erik Garrison
    Copyright © 2020-2024 Pjotr Prins

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


map<string, pair<vector<VariantAllele>,bool> > WfaVariant::wfa_parsedAlternates(
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
                /* DISABLED
                // flip the alt
                string alternate_rev = reverse_complement(alternate);
                // cout << "alt: " << alternate_rev << endl;
                alternate = alternate_rev;
                */
            }
            variantAlleles[alternate]; // create slot
            variantAlleles[alternate].second = is_inv;
            /*
            cerr << "comparing "
                 << " vs ref fwd " << rkmh::compare(alt_sketch, ref_sketch_fwd, invKmerLen)
                 << " vs rev " << rkmh::compare(alt_sketch, ref_sketch_rev, invKmerLen) << endl;
            */
        }

        if (is_inv) {
            // Do not realign and just annotate that there is an inversion
            variantAlleles[alternate].first.push_back(VariantAllele(ref, alternate, position));
        } else {
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
    }

    return variantAlleles;
}

/*
@@ Post-process alleles to reduce the set and normalise counts.
 */

#define ALLELE_NULL -1

void Variant::reduceAlleles(
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
            for (auto& va : v.alt) {
                va += a->ref.substr(v.ref.size());
            }
            v.ref = a->ref;
        }
        v.alt.push_back(a->alt);

        int alleleIndex = v.alt.size();
        for (const auto& originalIndex : originalIndexes) {
            unpackedAlleleIndexes[originalIndex][v.position] = alleleIndex;
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
        for (const auto& genotypeStr : genotypeStrs) {
            int i;
            if (!convert(genotypeStr, i)) {
                genotypeIndexes.push_back(ALLELE_NULL);
            } else {
                genotypeIndexes.push_back(i);
            }
        }
        map<long unsigned int, vector<int> > positionIndexes;
        for (const auto oldIndex : genotypeIndexes) {
            for (const auto& v : variants) {
                const long unsigned int& p = v.first;
                if (oldIndex == 0) { // reference
                    positionIndexes[p].push_back(0);
                } else {
                    positionIndexes[p].push_back(unpackedAlleleIndexes[oldIndex][p]);
                }
            }
        }
        for (auto& v : variants) {
            Variant& variant = v.second;
            vector<int>& gtints = positionIndexes[v.first];
            vector<string> gtstrs;
            for (const auto i : gtints) {
                if (i != ALLELE_NULL) {
                    gtstrs.push_back(convert(i));
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


} // namespace
