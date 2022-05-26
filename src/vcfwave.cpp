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
    -n, --nextgen           next gen mode.
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
    bool nextGen  = false;
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
                {"nextgen", no_argument, 0, 'n'},
                {"debug", no_argument, 0, 'd'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "ndhkt:L:p:t:K:I:f:",
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

        case 'n':
            nextGen = true;
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

        if (nextGen) {

            struct trackinfo {
                size_t pos0 = 0;
                string ref0, alt0, ref1, algn;
                size_t pos1 = 0;
                size_t altidx;
                int relpos;
                int AC,AF,AN;
                bool is_rev = false;
            };
            typedef map<string, trackinfo> TrackInfo;
            TrackInfo unique;

            // for (auto k: varAlleles) {
            for (const auto [alt0, wfvalue] : varAlleles) {
                // auto alt0 = k.first;
                // auto wfvalue = k.second;
                bool is_rev = wfvalue.second;
                for (auto wfmatch: wfvalue.first) {
                    auto ref = wfmatch.ref;
                    auto aligned = wfmatch.alt;
                    auto wfpos = wfmatch.position;
                    int alt_index,AC,AN = -1;
                    double AF = 0.0;
                    string wftag = alt0+":"+to_string(wfpos)+":"+ref+"/"+aligned;
                    if (var.ref == aligned) {
                        cout << "EQ: ";
                    }
                    else {
                        auto index = [](vector<string> v, string allele) {
                            auto it = find(v.begin(), v.end(), allele);
                            return (it == v.end() ? throw std::runtime_error("Unexpected value error for allele "+allele) : it - v.begin() );
                        };
                        alt_index = index(var.alt,alt0); // throws error if missing
                        AC = stoi(var.info["AC"].at(alt_index));
                        AF = stod(var.info["AF"].at(alt_index));
                        AN = stoi(var.info["AN"].at(0));
                    }
                    auto relpos = wfpos - var.position;
                    auto u = &unique[wftag];
                    u->pos0 = var.position;
                    u->ref0 = var.ref;
                    u->alt0 = alt0;
                    u->algn = aligned;
                    u->pos1 = wfpos;
                    u->altidx = alt_index;
                    u->relpos = relpos;
                    u->AC = AC;
                    u->AF = AF;
                    u->AN = AN;
                    u->is_rev = is_rev;
                    cout << wftag << " " << unique.size() << endl;
                }
            }
            cerr << "WIP " << endl;
        }
        else {
            var.legacy_reduceAlleles(
                varAlleles,
                variantFile,
                var,
                parseFlag,
                keepInfo,
                keepGeno,
                debug);
        }
    }

    return 0;
}
