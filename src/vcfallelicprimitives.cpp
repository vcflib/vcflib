/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2022 Erik Garrison
    Copyright © 2020-2022 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "bindings/cpp/WFAligner.hpp"
using namespace wfa;

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
    cerr << "usage: " << argv[0] << " [options] [file]" << endl
         << endl
         << "Realign reference and alternate alleles with WFA, parsing out the primitive alleles" << endl
         << "into multiple VCF records. New records have IDs that reference the source record ID." << endl
         << "Genotypes are handled. Deletion alleles will result in haploid (missing allele) genotypes." << endl
         << endl
         << "options:" << endl
         << "    -a, --algorithm TYPE    Choose algorithm (default) Wave front or (obsolete) Smith-Waterman" << endl
         << "                            [WF|SW] algorithm" << endl
         << "    -m, --use-mnps          Retain MNPs as separate events (default: false)." << endl
         << "    -t, --tag-parsed FLAG   Annotate decomposed records with the source record position" << endl
         << "                            (default: ORIGIN)." << endl
         << "    -L, --max-length LEN    Do not manipulate records in which either the ALT or" << endl
         << "                            REF is longer than LEN (default: unlimited)." << endl
         << "    -k, --keep-info         Maintain site and allele-level annotations when decomposing." << endl
         << "                            Note that in many cases, such as multisample VCFs, these won't" << endl
         << "                            be valid post-decomposition.  For biallelic loci in single-sample" << endl
         << "                            VCFs, they should be usable with caution." << endl
         << "    -d, --debug             debug mode." << endl;
    cerr << endl << "Type: transformation" << endl << endl;
    exit(0);
}

int main(int argc, char** argv) {

    bool includePreviousBaseForIndels = true;
    bool useMNPs = false;
    string parseFlag = "ORIGIN";
    string algorithm = "WF";
    int maxLength = 0;
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
           var.parsedAlternates(includePreviousBaseForIndels, useMNPs,
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
                if (debug) cerr << a.first << " " << va.repr << endl;
                alleles.insert(va); // only inserts first unique allele and ignores next ones
            }
        }

        int altcount = 0;
        for (auto a: alleles) {
            if (a.ref != a.alt) {
                ++altcount;
                if (debug) cerr << altcount << "$" << a.repr << endl;
            }
        }

        if (altcount == 1 && var.alt.size() == 1 && var.alt.front().size() == 1) { // if biallelic SNP
            cout << var << endl;
            continue;
        }

        // collect variant allele indexed membership
        map<string, vector<int> > variantAlleleIndexes; // from serialized VariantAllele to indexes
        for (auto a: varAlleles) {
            int index = var.altAlleleIndexes[a.first] + 1; // make non-relative
            for (auto va: a.second) {
                variantAlleleIndexes[va.repr].push_back(index);
            }
        }

        struct var_info_t {
            double freq = 0;
            int count = 0;
            map<string, string> info;
        };
        map<VariantAllele, var_info_t> alleleStuff;

        bool hasAf = false;
        if (var.info.find("AF") != var.info.end()) {
            hasAf = true;
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                vector<VariantAllele>& vars = varAlleles[*a];
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
            for (map<string, vector<string> >::iterator infoit = var.info.begin();
                 infoit != var.info.end(); ++infoit) {
                string key = infoit->first;
                for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                    vector<VariantAllele>& vars = varAlleles[*a];
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

        map<long unsigned int, Variant> variants;
        int varidx = 0;
        //vector<Variant> variants;
        for (set<VariantAllele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            if (a->ref == a->alt) {
                // ref allele
                continue;
            }
            vector<int>& originalIndexes = variantAlleleIndexes[a->repr];
            string type;
            int len = 0;
            if (a->ref.at(0) == a->alt.at(0)) { // well-behaved indels
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
            if (v.id.empty()) {
                v.id = var.id + "_" + std::to_string(++varidx);
            }
            //v.format = var.format;
            vector<string> gtonlyformat;
            gtonlyformat.push_back("GT");
            v.format = gtonlyformat;
            v.info["TYPE"].push_back(type);
            v.info["LEN"].push_back(convert(len));
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
            }
            // add null allele
            unpackedAlleleIndexes[ALLELE_NULL][v.position] = ALLELE_NULL;

        }

        // handle deletions
        for (set<VariantAllele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            int len = 0;
            if (a->ref.at(0) == a->alt.at(0)
                && a->ref.size() > a->alt.size()) {
                len = a->ref.size() - a->alt.size();
            } else {
                continue;
            }
            assert(len > 0);
            // nullify all the variants inside of the deletion range
            vector<int>& originalIndexes = variantAlleleIndexes[a->repr];
            auto begin = variants.upper_bound(a->position);
            auto end = variants.upper_bound(a->position + a->ref.size());
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

        //for (vector<Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
        for (map<long unsigned int, Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
            cout << v->second << endl;
        }
    }

    return 0;

}
