/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2022 Erik Garrison
    Copyright © 2020-2022 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include "convert.h"
#include "join.h"
#include "split.h"
#include <set>
#include <algorithm>
#include <getopt.h>

using namespace std;
using namespace vcflib;

#define ALLELE_NULL -1

double convertStrDbl(const string& s) {
    double r;
    convert(s, r);
    return r;
}

void printSummary(char** argv) {
    std::string text = R"(
usage: ./vcfallelicprimitives [options] [file]

WARNING: this tool is considered legacy and is only retained for older
workflows.  It will emit a warning!  Even though it can use the WFA
you should use [vcfwave](./vcfwave.md) instead.

Realign reference and alternate alleles with WFA or SW, parsing out
the primitive alleles into multiple VCF records. New records have IDs
that reference the source record ID.  Genotypes are handled. Deletion
alleles will result in haploid (missing allele) genotypes.

options:
    -a, --algorithm TYPE    Choose algorithm (default) Wave front or (obsolete)
                            Smith-Waterman [WF|SW] algorithm
    -m, --use-mnps          Retain MNPs as separate events (default: false).
    -t, --tag-parsed FLAG   Annotate decomposed records with the source record
                            position (default: ORIGIN).
    -L, --max-length LEN    Do not manipulate records in which either the ALT or
                            REF is longer than LEN (default: unlimited).
    -k, --keep-info         Maintain site and allele-level annotations when
                            decomposing.  Note that in many cases,
                            such as multisample VCFs, these won't be
                            valid post decomposition.  For biallelic
                            loci in single-sample VCFs, they should be
                            used with caution.
    -d, --debug             debug mode.

Type: transformation
)";

    cerr << text;
    exit(0);
}

int main(int argc, char** argv) {

    bool includePreviousBaseForIndels = true;
    bool useMNPs = false;
    string parseFlag = "ORIGIN";
    string algorithm = "WF";
    int maxLength = 0;
    size_t global_max_length = 0;
    bool keepInfo = false;
    bool keepGeno = false;
    bool useWaveFront = true;
    bool debug    = false;

    VariantCallFile variantFile;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"algorithm", required_argument, 0, 'a'},
                {"use-mnps", no_argument, 0, 'm'},
                {"max-length", required_argument, 0, 'L'},
                {"tag-parsed", required_argument, 0, 't'},
                {"keep-info", no_argument, 0, 'k'},
                {"keep-geno", no_argument, 0, 'g'},
                {"debug", no_argument, 0, 'd'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "dhmkt:L:a:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

        case 'a':
            algorithm = optarg;
            useWaveFront = (algorithm == "WF");
            break;


	    case 'm':
            useMNPs = true;
            break;

    case 'k':
            keepInfo = true;
            break;

	    case 'g':
            keepGeno = true;
            break;

	    case 'd':
            debug = true;
            break;

        case 'h':
            printSummary(argv);
            break;

	    case 't':
            parseFlag = optarg;
            break;

        case 'L':
            maxLength = atoi(optarg);
            break;

        case '?':
            printSummary(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (optind < argc) {
        string filename = argv[optind];
        variantFile.open(filename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    variantFile.addHeaderLine("##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">");
    variantFile.addHeaderLine("##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">");
    if (!parseFlag.empty()) {
      if (useWaveFront)
        variantFile.addHeaderLine("##INFO=<ID="+parseFlag+",Number=1,Type=String,Description=\"Decomposed from a complex record using vcflib vcfallelicprimitives and alignment with WFA2-lib.\">");
      else
        variantFile.addHeaderLine("##INFO=<ID="+parseFlag+",Number=1,Type=String,Description=\"Decomposed from a complex record using vcflib vcfallelicprimitives and alignment with obsolete SW.\">");
    }
    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        // we can't decompose *1* bp events, these are already in simplest-form whether SNPs or indels
        // we also don't handle anything larger than maxLength bp
        int max_allele_length = 0;
        for (auto allele: var.alt) {
          if (debug) cerr << allele << ":" << allele.length() << "," << max_allele_length << endl;
          global_max_length = max(allele.length(),global_max_length);
          if (allele.length() >= max_allele_length) {
             max_allele_length = allele.length();
             // cerr << max_allele_length << endl;
          }
        }

        if ((maxLength && max_allele_length > maxLength) || max_allele_length == 1 ||
            (var.alt.size() == 1 &&
             (var.ref.size() == 1 || (maxLength && var.ref.size() > maxLength)))) {
            // nothing to do
            cout << var << endl;
            continue;
        }

        // for each parsedAlternate allele, get the position
        // and build a new vcf record for that position
        // unless we are already at the position!
        // take everything which is unique to that allele (records) and append it to the new record
        // then handle genotypes; determine the mapping between allelic primitives and convert to phased haplotypes
        // this means taking all the parsedAlternates and, for each one, generating a pattern of allele indexes corresponding to it

        // this code does an O(n^2) alignment of the ALTs
        map<string, vector<VariantAllele> > varAlleles =
           var.legacy_parsedAlternates(includePreviousBaseForIndels, useMNPs,
                                false, // bool useEntropy = false,
                                10.0f, // float matchScore = 10.0f,
                                -9.0f, // float mismatchScore = -9.0f,
                                15.0f, // float gapOpenPenalty = 15.0f,
                                6.66f, // float gapExtendPenalty = 6.66f,
                                0.0f,  // float repeatGapExtendPenalty = 0.0f,
                                "",    // string flankingRefLeft = "",
                                "",    // string flankingRefRight = "",
                                useWaveFront,
                                debug);  // bool debug=false

        set<VariantAllele> alleles;
        // collect unique alleles
        for (auto a: varAlleles) {
            for (auto va: a.second) {
                if (debug) cerr << a.first << " " << va << endl;
                alleles.insert(va); // only inserts first unique allele and ignores if two are the same
            }
        }

        // count unique alleles
        int altcount = 0;
        for (auto a: alleles) {
            if (a.ref != a.alt) {
                ++altcount;
                if (debug) cerr << altcount << "$" << a << endl;
            }
        }

        if (altcount == 1 && var.alt.size() == 1 && var.alt.front().size() == 1) {
            // stop processing if biallelic SNP
            cout << var << endl;
            continue;
        }

        // collect variant allele indexed membership
        map<VariantAllele, vector<int> > variantAlleleIndexes; // from serialized VariantAllele to indexes
        for (auto a: varAlleles) {
            int index = var.altAlleleIndexes[a.first] + 1; // make non-relative
            for (auto va: a.second) {
                variantAlleleIndexes[va].push_back(index);
            }
        }

        // Here we correct the fields in the VCF record
        struct var_info_t {
            double freq = 0;
            int count = 0;
            map<string, string> info;
        };
        map<VariantAllele, var_info_t> alleleStuff;

        // AF: Allele frequency field
        bool hasAf = false;
        if (var.info.find("AF") != var.info.end()) {
            hasAf = true;
            for (auto a: var.alt) {
                auto& vars = varAlleles[a];
                for (auto va: vars) {
                    double freq;
                    try {
                        convert(var.info["AF"].at(var.altAlleleIndexes[a]), freq);
                        alleleStuff[va].freq += freq;
                    } catch (...) {
                        cerr << "vcfallelicprimitives WARNING: AF does not have count == alts @ "
                             << var.sequenceName << ":" << var.position << endl;
                    }
                }
            }
        }

        // AC: Allele count field
        bool hasAc = false;
        if (var.info.find("AC") != var.info.end()) {
            hasAc = true;
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                vector<VariantAllele>& vars = varAlleles[*a];
                for (vector<VariantAllele>::iterator va = vars.begin(); va != vars.end(); ++va) {
                    int freq;
                    try {
                        convert(var.info["AC"].at(var.altAlleleIndexes[*a]), freq);
                        alleleStuff[*va].count += freq;
                    } catch (...) {
                        cerr << "vcfallelicprimitives WARNING: AC does not have count == alts @ "
                             << var.sequenceName << ":" << var.position << endl;
                    }
                }
            }
        }

        if (keepInfo) {
            for (auto infoit: var.info) {
                string key = infoit.first;
                for (auto a: var.alt) {
                    vector<VariantAllele>& vars = varAlleles[a];
                    for (auto va: vars) {
                        string val;
                        vector<string>& vals = var.info[key];
                        if (vals.size() == var.alt.size()) { // allele count for info
                            val = vals.at(var.altAlleleIndexes[a]);
                        } else if (vals.size() == 1) { // site-wise count
                            val = vals.front();
                        } // don't handle other multiples... how would we do this without going crazy?
                        if (!val.empty()) {
                            alleleStuff[va].info[key] = val;
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

        map<long unsigned int, Variant> variants;
        int varidx = 0;
        for (auto a: alleles) {
            const auto ref = a.ref;
            const auto alt = a.alt;

            if (a.ref == a.alt)
                continue;

            vector<int>& originalIndexes = variantAlleleIndexes[a];
            string type;
            int len = 0;
            if (ref.at(0) == alt.at(0)) { // well-behaved indels
                if (ref.size() > alt.size()) {
                    type = "del";
                    len = ref.size() - alt.size();
                    // special case
                    // a deletion implies we should be ALLELE_NULL on this haplotype
                    // until the end of the deletion
                    // save the range in a new map which we'll iterate over
                    for (auto i : originalIndexes) {
                        // TODO check if it should be len
                        //auto d = (*deletions)[i];
                        //d.push_back(make_pair(0, 0));
                    }
                } else if (ref.size() < alt.size()) {
                    len = alt.size() - ref.size();
                    type = "ins";
                }
            } else {
                if (ref.size() == alt.size()) {
                    len = ref.size();
                    if (ref.size() == 1) {
                        type = "snp";
                    } else {
                        type = "mnp";
                    }
                } else {
                    len = abs((int) ref.size() - (int) alt.size());
                    type = "complex";
                }
            }

            if (variants.find(a.position) == variants.end()) {
                Variant newvar(variantFile);
                variants[a.position] = newvar;
            }

            Variant& v = variants[a.position]; // guaranteed to exist

            if (!parseFlag.empty()) {
                v.info[parseFlag].push_back(var.sequenceName + ":" + std::to_string(var.position));
            }
            v.quality = var.quality;
            v.filter = var.filter;
            v.infoOrderedKeys= var.infoOrderedKeys;
            if (v.id.empty()) {
                v.id = var.id + "_" + std::to_string(++varidx);
            }
            //v.format = var.format;
            vector<string> gtonlyformat;
            gtonlyformat.push_back("GT");
            v.format = gtonlyformat;
            v.infoOrderedKeys.push_back("LEN");
            v.infoOrderedKeys.push_back("ORIGIN");
            v.infoOrderedKeys.push_back("TYPE");
            v.info["TYPE"].push_back(type);
            v.info["LEN"].push_back(convert(len));
            if (hasAf) {
                v.info["AF"].push_back(convert(alleleStuff[a].freq));
            }
            if (hasAc) {
                v.info["AC"].push_back(convert(alleleStuff[a].count));
            }
            if (keepInfo) {
                for (auto infoit: var.info) {
                    string key = infoit.first;
                    if (key != "AF" && key != "AC" && key != "TYPE" && key != "LEN") { // don't clobber previous
                        v.info[key].push_back(alleleStuff[a].info[key]);
                    }
                }
            }

            // now, keep all the other infos if we are asked to
            v.sequenceName = var.sequenceName;
            v.position = a.position; // ... by definition, this should be == if the variant was found
            if (v.ref.size() < ref.size()) {
                for (auto &va: v.alt) {
                    va += ref.substr(v.ref.size());
                }
                v.ref = ref;
            }
            v.alt.push_back(alt);

            int alleleIndex = v.alt.size();
            for (auto i: originalIndexes) {
                unpackedAlleleIndexes[i][v.position] = alleleIndex;
            }
            // add null allele
            unpackedAlleleIndexes[ALLELE_NULL][v.position] = ALLELE_NULL;

        }

        // handle deletions
        for (auto a: alleles) {
            const auto ref = a.ref;
            const auto alt = a.alt;
            int len = 0;
            if (ref.at(0) == alt.at(0)
                && ref.size() > alt.size()) {
                len = ref.size() - alt.size();
            } else {
                continue;
            }
            assert(len > 0);
            // nullify all the variants inside of the deletion range
            vector<int>& originalIndexes = variantAlleleIndexes[a];
            auto begin = variants.upper_bound(a.position);
            auto end = variants.upper_bound(a.position + ref.size());
            for (auto i : originalIndexes) {
                for (auto x = begin; x != end; ++x) {
                    unpackedAlleleIndexes[i][x->second.position] = ALLELE_NULL;
                }
            }
        }

        // genotypes
        for (auto s: var.sampleNames) {
            string& sampleName = s;
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
            for (auto gs: genotypeStrs) {
                int i;
                if (!convert(gs, i)) {
                    genotypeIndexes.push_back(ALLELE_NULL);
                } else {
                    genotypeIndexes.push_back(i);
                }
            }
            map<long unsigned int, vector<int> > positionIndexes;
            for (auto g: genotypeIndexes) {
                int oldIndex = g;
                for (map<long unsigned int, Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
                    const long unsigned int& p = v->first;
                    if (oldIndex == 0) { // reference
                        positionIndexes[p].push_back(0);
                    } else {
                        positionIndexes[p].push_back(unpackedAlleleIndexes[oldIndex][p]);
                    }
                }
            }
            for (auto &v: variants) {
                Variant& variant = v.second;
                vector<int>& gtints = positionIndexes[v.first];
                vector<string> gtstrs;
                for (auto i: gtints) {
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
    if (global_max_length > 100) {
        cerr << "WARNING: this tool is considered legacy. This VCF file contains sequences longer than 100 bps. Please use vcfwave instead!" << endl;
    }

    return 0;

}
