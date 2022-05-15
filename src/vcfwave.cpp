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
#include <omp.h>
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
usage: vcfwave [options] [file]

Realign reference and alternate alleles with WFA, parsing out the primitive alleles
into multiple VCF records. New records have IDs that reference the source record ID.
Genotypes are handled. Deletions generate haploid/missing genotypes at overlapping sites.

options:
    -p, --wf-params PARAMS  use the given BiWFA params (default: 0,19,39,3,81,1)
                            format=match,mismatch,gap1-open,gap1-ext,gap2-open,gap2-ext
    -f, --tag-parsed FLAG   Annotate decomposed records with the source record position
                            (default: ORIGIN).
    -L, --max-length LEN    Do not manipulate records in which either the ALT or
                            REF is longer than LEN (default: unlimited).
    -K, --inv-kmer K        Length of k-mer to use for inversion detection sketching (default: 17).
    -I, --inv-min LEN       Minimum allele length to consider for inverted alignment (default: 64).
    -k, --keep-info         Maintain site and allele-level annotations when decomposing.
                            Note that in many cases, such as multisample VCFs, these won't
                            be valid post-decomposition.  For biallelic loci in single-sample
                            VCFs, they should be usable with caution.
    -t, --threads N         use this many threads for variant decomposition
    -d, --debug             debug mode.

Type: transformation
)";

    cerr << text;

    exit(0);
}

int main(int argc, char** argv) {

    bool includePreviousBaseForIndels = true;
    bool useMNPs = true;
    string parseFlag = "ORIGIN";
    string algorithm = "WF";
    string paramString = "0,19,39,3,81,1";
    int maxLength = 0;
    bool keepInfo = false;
    bool keepGeno = false;
    bool useWaveFront = true;
    bool debug    = false;
    int thread_count = 1;
    int inv_sketch_kmer = 17;
    int min_inv_len = 64;

    VariantCallFile variantFile;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"wf-params", required_argument, 0, 'p'},
                {"max-length", required_argument, 0, 'L'},
                {"inv-kmer", required_argument, 0, 'K'},
                {"inv-min", required_argument, 0, 'I'},
                {"tag-parsed", required_argument, 0, 'f'},
                {"keep-info", no_argument, 0, 'k'},
                {"keep-geno", no_argument, 0, 'g'},
                {"threads", required_argument, 0, 't'},
                {"debug", no_argument, 0, 'd'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "dhkt:L:p:t:K:I:f:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

	    case 'k':
            keepInfo = true;
            break;

	    case 'g':
            keepGeno = true;
            break;

        case 'p':
            paramString = optarg;
            break;

        case 't':
            thread_count = atoi(optarg);
            break;

	    case 'd':
            debug = true;
            break;

        case 'h':
            printSummary(argv);
            break;

	    case 'f':
            parseFlag = optarg;
            break;

        case 'L':
            maxLength = atoi(optarg);
            break;

        case 'K':
            inv_sketch_kmer = atoi(optarg);
            break;

        case 'I':
            min_inv_len = atoi(optarg);
            break;

        case '?':
            printSummary(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    omp_set_num_threads(thread_count);

    if (optind < argc) {
        string filename = argv[optind];
        variantFile.open(filename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    // parse the alignment parameters
    vector<string> p_str = split(paramString, ',');
    vector<int> p;
    for (auto& s : p_str) { p.push_back(atoi(s.c_str())); }

    auto wfa_params = wavefront_aligner_attr_default;
    wfa_params.memory_mode = wavefront_memory_ultralow; // note this is overridden in Variant.cpp
    wfa_params.distance_metric = gap_affine_2p;
    wfa_params.affine2p_penalties.match = p[0];
    wfa_params.affine2p_penalties.mismatch = p[1];
    wfa_params.affine2p_penalties.gap_opening1 = p[2];
    wfa_params.affine2p_penalties.gap_extension1 = p[3];
    wfa_params.affine2p_penalties.gap_opening2 = p[4];
    wfa_params.affine2p_penalties.gap_extension2 = p[5];
    wfa_params.alignment_scope = compute_alignment;

    variantFile.addHeaderLine("##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">");
    variantFile.addHeaderLine("##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">");
    variantFile.addHeaderLine("##INFO=<ID="+parseFlag+",Number=1,Type=String,Description=\"Decomposed from a complex record using vcflib vcfwave and alignment with WFA2-lib.\">");
    variantFile.addHeaderLine("##INFO=<ID=INV,Number=A,Type=String,Description=\"Count of haplotypes which are aligned in the inverted orientation using vcflib vcfwave.\">");
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

        map<string, pair<vector<VariantAllele>, bool> > varAlleles =
           var.parsedAlternates(includePreviousBaseForIndels, useMNPs,
                                false, // bool useEntropy = false,
                                "",    // string flankingRefLeft = "",
                                "",    // string flankingRefRight = "",
                                &wfa_params,
                                inv_sketch_kmer,
                                min_inv_len,
                                debug);  // bool debug=false

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
            continue;
        }

        // collect variant allele indexed membership
        map<VariantAllele, vector<int> > variantAlleleIndexes; // from serialized VariantAllele to indexes
        for (auto a: varAlleles) {
            int index = var.altAlleleIndexes[a.first] + 1; // make non-relative
            for (auto va: a.second.first) {
                variantAlleleIndexes[va].push_back(index);
            }
        }

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
        //vector<Variant> variants;
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

        // handle deletions
        for (set<VariantAllele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            int len = 0;
            if (a->ref.size() && a->alt.size() && a->ref.at(0) == a->alt.at(0)
                && a->ref.size() > a->alt.size()) {
                len = a->ref.size() - a->alt.size();
            } else {
                continue;
            }
            assert(len > 0);
            // nullify all the variants inside of the deletion range
            vector<int>& originalIndexes = variantAlleleIndexes[*a];
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
