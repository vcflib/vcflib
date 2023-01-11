/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2023 Erik Garrison
    Copyright © 2020-2023 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include <algorithm>
#include <string>
#include <set>
#include <utility>
#include <vector>

#include <getopt.h>
#include <omp.h>

#include "Variant.h"
#include "convert.h"
#include "join.h"
#include "split.h"
#include "progress.h"

using namespace std;
using namespace vcflib;

#define ALLELE_NULL -1
#define ALLELE_NULL2 -200 // large number brings out issues

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
    --quiet                 no progress bar
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
    bool nextGen  = true;
    bool quiet    = false;
    bool debug    = false;

    int thread_count = 1;
    int inv_sketch_kmer = 17;
    int min_inv_len = 64;

    VariantCallFile variantFile;

    // int y = MIN(2,3);
    // auto a = DNA_CHAR_A;
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
                {"quiet", no_argument, 0, 'q'},
                {"debug", no_argument, 0, 'd'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "nqdhkt:L:p:t:K:I:f:",
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

        case 'q':
            quiet = true;
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

    off_t file_size = -1;

    if (optind < argc) {
        string filename = argv[optind];
        variantFile.open(filename);
        file_size = get_file_size(filename.c_str());
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
    double amount = 0.0;
    uint64_t start = get_timestamp();

    if (!quiet)
        cerr << "vcfwave processing VCF file..." << endl;
    while (variantFile.getNextVariant(var)) {

        amount = (double)variantFile.file_pos()/(double)file_size;
        // cerr << file_size << "," << variantFile.file_pos() << "=" << amount << endl;
        if (!quiet && variantFile.file_pos() >= 0 && file_size >= 0)
            print_progress(amount*100, start);

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
            // The following section post-process the results of wavefront and
            // updates AC, AF and genotype values
            typedef vector<int> Genotypes;
            typedef vector<Genotypes> RecGenotypes;
            struct trackinfo {
                size_t pos0 = 0;
                string ref0, alt0, ref1, algn;
                size_t pos1 = 0;
                size_t altidx;
                int relpos;
                int AC=0,AN=0;
                double AF=0.0;
                string AT;
                int size = -99;
                bool is_inv = false;
                string type;
                string origin;
                RecGenotypes genotypes;
            };
            typedef map<string, trackinfo> TrackInfo;
            TrackInfo unique; // Track all alleles

            // Unpack wavefront results and set values for each unique allele
            for (const auto [alt0, wfvalue] : varAlleles) {
                bool is_inv = wfvalue.second;
                for (auto wfmatch: wfvalue.first) {
                    auto ref = wfmatch.ref;
                    auto aligned = wfmatch.alt;
                    auto wfpos = wfmatch.position;
                    int alt_index,AC,AN = -1;
                    string AT;
                    double AF = 0.0;
                    string wftag = alt0+":"+to_string(wfpos)+":"+ref+"/"+aligned;
                    if (var.ref != aligned) {
                        auto index = [&](vector<string> v, string allele) {
                            auto check = (is_inv ? reverse_complement(allele) : allele);
                            auto it = find(v.begin(), v.end(), check);
                            return (it == v.end() ? throw std::runtime_error("Unexpected value error for allele (inv="+to_string(is_inv)+ " " +check) : it - v.begin() );
                        };
                        alt_index = index(var.alt,alt0); // throws error if missing
                        AC = stoi(var.info["AC"].at(alt_index));
                        AF = stod(var.info["AF"].at(alt_index));
                        AT = var.info["AT"].at(alt_index);
                        AN = stoi(var.info["AN"].at(0));
                    }
                    auto relpos = wfpos - var.position;
                    auto u = &unique[wftag];
                    u->pos0 = var.position;
                    u->ref0 = var.ref;
                    u->alt0 = alt0;
                    u->ref1 = ref;
                    u->algn = aligned;
                    u->pos1 = wfpos;
                    u->altidx = alt_index;
                    u->relpos = relpos;
                    u->AC = AC;
                    u->AF = AF;
                    u->AN = AN;
                    u->AT = AT;
                    u->is_inv = is_inv;
                }
            }
            // Collect genotypes for every allele from the main record. This code is
            // effectively mirrored in Python in realign.py:
            RecGenotypes genotypes;
            auto samples = var.samples;
            for (auto sname: var.sampleNames) {
                auto genotype1 = samples[sname]["GT"].front();
                vector<string> genotypeStrs = split(genotype1, "|/");
                Genotypes gts;
                std::transform(genotypeStrs.begin(), genotypeStrs.end(), std::back_inserter(gts), [](auto n){ return (n == "." ? ALLELE_NULL2 : stoi(n)); });
                genotypes.push_back(gts);
            }
            // Now plug in the new indices for listed genotypes
            for (auto [tag,aln]: unique) {
                RecGenotypes aln_genotypes = genotypes; // make a copy
                auto altidx1 = aln.altidx+1;
                for (auto &gt: aln_genotypes) {
                    int i = 0;
                    for (auto g: gt) {
                        if (g == altidx1)
                            gt[i] = 1; // one genotype in play
                        else
                            if (g != ALLELE_NULL2) gt[i] = 0;
                        i++;
                    }
                }
                unique[tag].genotypes = aln_genotypes;
            }

            // Merge records that describe the exact same variant (in
            // earlier jargon a 'primitive allele' in a new dict named
            // variants and adjust AC, AF and genotypes:
            TrackInfo track_variants;
            for (auto [key,v] : unique) {
                auto ref = v.ref1;
                auto aligned = v.algn;
                if (ref != aligned) {
                    auto ntag = to_string(v.pos1) + ":" + ref + "/" + aligned + "_" + to_string(v.is_inv);
                    if (track_variants.count(ntag)>0) { // this variant already exists
                        track_variants[ntag].AC += v.AC;
                        // Check AN number is equal so we can compute AF by addition
                        assert(track_variants[ntag].AN == v.AN);
                        track_variants[ntag].AF += v.AF;
                        // Merge genotypes if they come from different alleles
                        if (v.altidx != track_variants[ntag].altidx) {
                            auto track_genotypes = track_variants[ntag].genotypes;
                            for (auto sample: v.genotypes) { // all samples
                                int i = 0;
                                for (auto g: sample) { // all genotypes
                                    if (g && sample[i]>0)
                                        sample[i] = 1; // always one genotype in play
                                    i++;
                                }
                            }
                        }
                    }
                    else {
                        track_variants[ntag] = v;
                    }
                }
            }
            unique.clear();
            // The following section updates the INFO TYPE and INV field:
            // Adjust TYPE field to set snp/mnp/ins/del
            for (auto [key,v] : track_variants) {
                auto ref_len = v.ref1.length();
                auto aln_len = v.algn.length();
                string type;
                auto size = -99;
                if (aln_len < ref_len) {
                    type = "del";
                    size = ref_len - aln_len;
                }
                else if (aln_len > ref_len) {
                    type = "ins";
                    size = aln_len - ref_len;
                }
                else if (aln_len == ref_len) {
                    if (ref_len == 1)
                        type = "snp";
                    else
                        type = "mnp";
                    size = aln_len;
                }

                v.type = type;
                v.size = size;
                v.origin = var.sequenceName+":"+to_string(var.position);
                track_variants[key] = v;
            }
            // Here we correct for deletions - overlapping cals for SNP and MNP get nullified.
            for (auto [key,v]: track_variants) {
                if (v.type == "del") {
                    auto del_ref_len = v.ref1.length();
                    auto del_aln_len = v.algn.length();
                    auto del_pos1 = v.pos1;
                    auto del_size = v.size;
                    auto del_start_pos = del_pos1 + del_aln_len;
                    // Make a range from the start of the deletion to the end
                    auto check_range = make_tuple(del_start_pos, del_start_pos + del_size);
                    auto check_samples = v.genotypes;
                    for (auto [key2,var2]: track_variants) {
                        if (var2.type == "snp" || var2.type == "mnp") {
                            // for alignment check all SNPs/MNPs
                            auto pos1 = var2.pos1;
                            auto pos2 = pos1 + var2.size;
                            auto overlap = [] (unsigned int pos,tuple<unsigned int, unsigned int> range) {
                                auto start = get<0>(range);
                                auto end = get<1>(range);
                                return (pos >= start || pos <= end);
                            };
                            if (overlap(pos1,check_range) || overlap(pos2,check_range)) {
                                int i = 0;
                                for (auto &sample2: var2.genotypes) {
                                    auto del_sample = check_samples[i];
                                    auto find_del = find(del_sample.begin(), del_sample.end(), 1);
                                    bool nullify = !(find_del == del_sample.end());
                                    if (nullify) {
                                        for (auto &item: sample2) {
                                            item = ALLELE_NULL2;
                                        }
                                    }
                                    i++;
                                }

                            }
                            track_variants[key2] = var2;
                        }
                    }
                }
            }
            // The following section outputs all tracked alleles one by one:
            int ct = 0;
            for (auto [key,v]: track_variants) {
                ct++;
                Variant newvar(variantFile);
                newvar.sequenceName = var.sequenceName;
                newvar.position = v.pos1;
                newvar.id = var.id + "_" + to_string(ct);
                newvar.ref = v.ref1;
                newvar.alt.push_back(v.algn);
                newvar.quality = var.quality;
                newvar.info = var.info;
                newvar.infoOrderedKeys = var.infoOrderedKeys;

                vector<string> AT{ v.AT };
                vector<string> ORIGIN{ v.origin };
                vector<string> TYPE{ v.type };
                newvar.info["AC"] = vector<string>{ to_string(v.AC) };
                newvar.info["AF"] = vector<string>{ to_string(v.AF) };
                newvar.info["AN"] = vector<string>{ to_string(v.AN) };
                newvar.info["AT"] = AT;
                newvar.info[parseFlag] = ORIGIN;
                newvar.info["TYPE"] = TYPE;
                newvar.info["LEN"] = vector<string>{to_string(v.size)};
                newvar.info["INV"] = vector<string>{to_string(v.is_inv)};
                // set the output order of the new INFO fields:
                newvar.infoOrderedKeys.push_back("ORIGIN");
                newvar.infoOrderedKeys.push_back("LEN");
                newvar.infoOrderedKeys.push_back("INV");
                newvar.infoOrderedKeys.push_back("TYPE");
                // newvar.format = var.format;
                // newvar.sampleNames = var.sampleNames;
                // newvar.outputSampleNames = var.outputSampleNames;
                // newvar.samples = v.genotypeStrs;

                // Instead of using above format output we now simply print genotypes
                cout.precision(2);
                cout << newvar;
                cout << "\tGT";
                for (auto gts: v.genotypes) {
                    cout << "\t";
                    int idx = 0;
                    for (auto gt : gts) {
                        cout << (gt == ALLELE_NULL2 ? "." : to_string(gt));
                        if (idx < gts.size()-1) cout << "|";
                        idx++;
                    }
                }
                cout << endl;
            }
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

    if (!quiet) cerr << endl;

    return 0;
}
