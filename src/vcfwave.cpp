/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2024 Erik Garrison
    Copyright © 2020-2024 Pjotr Prins
    Copyright © 2024 Andrea Guarracino

    This software is published under the MIT License. See the LICENSE file.
*/

#include <algorithm>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <sstream>

#include <getopt.h>
#ifdef WFA_PARALLEL
#include <omp.h>
#endif

#include "Variant.h"
#include "vcf-wfa.h"
#include "convert.h"
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
    const std::string text = R"(
usage: vcfwave [options] [file]

Realign reference and alternate alleles with WFA, parsing out the
'primitive' alleles into multiple VCF records. New records have IDs that
reference the source record ID.  Genotypes/samples are handled
correctly. Deletions generate haploid/missing genotypes at overlapping
sites.

options:
    -p, --wf-params PARAMS  use the given BiWFA params (default: 0,19,39,3,81,1)
                            format=match,mismatch,gap1-open,gap1-ext,gap2-open,gap2-ext
    -f, --tag-parsed FLAG   Annotate decomposed records with the source record position
                            (default: ORIGIN).
    -L, --max-length LEN    Do not manipulate records in which either the ALT or
                            REF is longer than LEN (default: unlimited).
    -K, --inv-kmer K        Length of k-mer to use for inversion detection sketching (default: 17).
    -I, --inv-min LEN       Minimum allele length to consider for inverted alignment (default: 64).
    -t, --threads N         Decompose up to N records in parallel (default 1); output order
                            is preserved. Peak memory grows ~0.7GB per concurrent big site.
    -w, --wfa-threads M     Legacy: M cores inside each WFA alignment, one record at a time
                            (default 1). Mutually exclusive with -t.
    --quiet                 Do not display progress bar.
    -d, --debug             Debug mode.

Note the -k,--keep-info switch is no longer in use and ignored.

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
    bool nextGen  = true; // postprocessing is the default now
    bool quiet    = false;
    bool debug    = false;

    int thread_count = 1;
    int wfa_thread_count = 1;
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
                {"help", no_argument, nullptr, 'h'},
                {"wf-params", required_argument, nullptr, 'p'},
                {"max-length", required_argument, nullptr, 'L'},
                {"inv-kmer", required_argument, nullptr, 'K'},
                {"inv-min", required_argument, nullptr, 'I'},
                {"tag-parsed", required_argument, nullptr, 'f'},
                {"keep-info", no_argument, nullptr, 'k'},
                {"keep-geno", no_argument, nullptr, 'g'},
                {"threads", required_argument, nullptr, 't'},
                {"wfa-threads", required_argument, nullptr, 'w'},
                {"nextgen", no_argument, nullptr, 'n'},
                {"quiet", no_argument, nullptr, 'q'},
                {"debug", no_argument, nullptr, 'd'},
                {nullptr, 0, nullptr, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "nqdhkt:w:L:p:K:I:f:",
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

        case 'w':
            wfa_thread_count = atoi(optarg);
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

    // -t (parallel records, WFA single-threaded) and -w (one record at a time,
    // WFA multi-threaded) are mutually exclusive: combining them would require
    // nested OpenMP, and the WFA-threaded path never outperforms record-level
    // parallelism anyway. Keeping them exclusive means the -w path runs the
    // record loop single-threaded, so WFA's threads are not nested.
    if (thread_count > 1 && wfa_thread_count > 1) {
        cerr << "vcfwave: -t/--threads and -w/--wfa-threads are mutually exclusive "
             << "(use -t to parallelise across records, or -w for the legacy "
             << "per-alignment threading)." << endl;
        exit(1);
    }

    #ifdef WFA_PARALLEL
    omp_set_num_threads(thread_count);
    #endif

    off_t file_size = -1;

    if (optind < argc) {
        const string filename = argv[optind];
        variantFile.open(filename);
        file_size = get_file_size(filename.c_str());
    }
    else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    // parse the alignment parameters
    vector<string> p_str = split(paramString, ',');
    vector<int> p;
    p.reserve(p_str.size());
    for (const auto& s : p_str) { p.push_back(atoi(s.c_str())); }

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
    variantFile.addHeaderLine("##INFO=<ID=INV,Number=0,Type=Flag,Description=\"Inversion detected\">");
    cout << variantFile.header << endl;

    double amount = 0.0, prev_amount = 0.0;
    uint64_t start = get_timestamp();

    if (!quiet)
        cerr << "vcfwave " << VCFLIB_VERSION << " processing..." << endl;

    // Decompose one record. Returns true and writes the decomposed VCF output to
    // `out` for records that are actually realigned; returns false WITHOUT
    // writing for pass-through records (the driver renders those serially, so the
    // historical carry-over of cout's float precision is preserved byte-for-byte).
    // All configuration is captured by reference and only read; `var` and `out`
    // are per-call, so this is safe to run concurrently on distinct records.
    auto process_record = [&](WfaVariant& var, ostream& out) -> bool {

        // we can't decompose *1* bp events, these are already in simplest-form whether SNPs or indels
        // we also don't handle anything larger than maxLength bp
        int max_allele_length = 0;
        for (const auto& allele: var.alt) {
          if (debug) cerr << allele << ":" << allele.length() << "," << max_allele_length << endl;
          if (allele.length() >= max_allele_length) {
             max_allele_length = allele.length();
             // cerr << max_allele_length << endl;
          }
        }

        if ((maxLength && max_allele_length > maxLength) || max_allele_length == 1 ||
            (var.alt.size() == 1 &&
             (var.ref.size() == 1 || (maxLength && var.ref.size() > maxLength)))) {
            // nothing to do -- pass-through; rendered serially by the driver
            return false;
        }

        map<string, pair<vector<VariantAllele>, bool> > varAlleles =
           var.wfa_parsedAlternates(includePreviousBaseForIndels, useMNPs,
                                false, // bool useEntropy = false,
                                "",    // string flankingRefLeft = "",
                                "",    // string flankingRefRight = "",
                                &wfa_params,
                                inv_sketch_kmer,
                                min_inv_len,
                                wfa_thread_count, // cores per WFA alignment (-w); parallelism is mainly across records
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
                int AC=-1,AN=-1;
                double AF=-1;
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
            for (const auto& [alt0, wfvalue] : varAlleles) {               
                bool is_inv = wfvalue.second;
                for (const auto& wfmatch: wfvalue.first) {
                    const auto& ref = wfmatch.ref;
                    const auto& aligned = wfmatch.alt;
                    const auto wfpos = wfmatch.position;
                    int alt_index=-1,AC=-1,AN = -1;
                    string AT;
                    double AF = -1;
                    if (ref != aligned) {
                        auto index = [&](const vector<string>& v, const string& allele) {
                            //auto check = (is_inv ? reverse_complement(allele) : allele); DISABLED
                            const auto& check = allele;
                            auto it = find(v.begin(), v.end(), check);
                            return it == v.end() ? throw std::runtime_error("Unexpected value error for allele (inv="+to_string(is_inv)+ " " +check + ")") : it - v.begin();
                        };
                        alt_index = index(var.alt,alt0); // throws error if missing
                        if (var.info["AC"].size() > alt_index) {
                            AC = stoi(var.info["AC"].at(alt_index));
                        }
                        if (var.info["AF"].size() > alt_index) {
                            AF = stod(var.info["AF"].at(alt_index));
                        }
                        if (var.info["AT"].size() > alt_index) {
                            AT = var.info["AT"].at(alt_index);
                        }
                        if (var.info["AN"].size() > alt_index) {
                            AN = stoi(var.info["AN"].at(alt_index));
                        }
                    }
                    const auto relpos = wfpos - var.position;

                	const string wftag = alt0 + ":" + to_string(wfpos) + ":" + ref + "/" + aligned;
                    auto& u = unique[wftag];
                    u.pos0 = var.position;
                    u.ref0 = var.ref;
                    u.alt0 = alt0;
                    u.ref1 = ref;
                    u.algn = aligned;
                    u.pos1 = wfpos;
                    u.altidx = alt_index;
                    u.relpos = relpos;
                    u.AC = AC;
                    u.AF = AF;
                    u.AN = AN;
                    u.AT = AT;
                    u.is_inv = is_inv;
                }
            }
            // Collect genotypes for every allele from the main record. This code is
            // effectively mirrored in Python in realign.py:
            RecGenotypes genotypes;
            auto samples = var.samples;
            for (const auto& sname: var.sampleNames) {
                const auto& genotype1 = samples[sname]["GT"].front();
                vector<string> genotypeStrs = split(genotype1, "|/");
                Genotypes gts;
                std::transform(genotypeStrs.begin(), genotypeStrs.end(), std::back_inserter(gts), [](const auto& n){ return (n == "." ? ALLELE_NULL2 : stoi(n)); });
                genotypes.push_back(gts);
            }
            // Now plug in the new indices for listed genotypes
            for (auto& [_,aln]: unique) {
                RecGenotypes aln_genotypes = genotypes; // make a copy
                const auto altidx1 = aln.altidx+1;
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
                aln.genotypes = aln_genotypes;
            }

            // Merge records that describe the exact same variant (in
            // earlier jargon a 'primitive allele' in a new dict named
            // variants and adjust AC, AF and genotypes:
            TrackInfo track_variants;
            for (const auto& [_,v] : unique) {
                const auto& ref = v.ref1;
                const auto& aligned = v.algn;
                if (ref != aligned) {
                    auto ntag = to_string(v.pos1) + ":" + ref + "/" + aligned + "_" + to_string(v.is_inv) + "_"+ v.AT; 
                    if (track_variants.count(ntag)>0 && track_variants[ntag].AN == v.AN) { // this variant already exists
                        track_variants[ntag].AC += v.AC;
                        // Check AN number is equal so we can compute AF by addition
                        //assert(track_variants[ntag].AN == v.AN);
                        track_variants[ntag].AF += v.AF;
                    }
                    else {
                        track_variants[ntag] = v;
                    }
                }
            }
            unique.clear();
            // The following section updates the INFO TYPE and INV field:
            // Adjust TYPE field to set snp/mnp/ins/del
            for (auto& [_,v] : track_variants) {
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
            }
            // Here we correct for deletions - overlapping cals for SNP and MNP get nullified.
            for (const auto& [_,v]: track_variants) {
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
                                return (pos >= start && pos < end);
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
            for (const auto& [_,v]: track_variants) {
                ct++;
                Variant newvar(variantFile);
                newvar.sequenceName = var.sequenceName;
                newvar.position = v.pos1;

                // inversions are not realigned anymore, so they can't generate multiple entries
                newvar.id = v.is_inv ? var.id : var.id + "_" + to_string(ct);
                
                newvar.ref = v.ref1;
                newvar.alt.push_back(v.algn);
                newvar.quality = var.quality;
                newvar.info = var.info;
                newvar.infoOrderedKeys = var.infoOrderedKeys;

                vector<string> AT{ v.AT };
                vector<string> TYPE{ v.type };
                if (v.AC > -1) {
                    newvar.info["AC"] = vector<string>{ to_string(v.AC) };
                }
                if (v.AF > -1) {
                    newvar.info["AF"] = vector<string>{ to_string(v.AF) };
                }
                if (v.AN > -1) {
                    newvar.info["AN"] = vector<string>{ to_string(v.AN) };
                }
                if (v.AT.find_first_not_of(' ') != std::string::npos) {
                    newvar.info["AT"] = AT; // there is a non-space character
                }

                // Inversions are not decomposed anymore, so there is no need to specify the ORIGIN
                if (!v.is_inv) {
                    vector<string> ORIGIN{ v.origin };
                    newvar.info[parseFlag] = ORIGIN;
                }
                newvar.info["TYPE"] = TYPE;
                newvar.info["LEN"] = vector<string>{to_string(v.size)};

                // Emit INV=YES if the variant is an inversion
                if (v.is_inv) {
                    newvar.info["INV"] = vector<string>{"YES"};
                }

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
                out.precision(2);
                out << newvar;
                out << "\tGT";
                for (const auto& gts: v.genotypes) {
                    out << "\t";
                    int idx = 0;
                    for (auto gt : gts) {
                        out << (gt == ALLELE_NULL2 ? "." : to_string(gt));
                        if (idx < gts.size()-1) out << "|";
                        idx++;
                    }
                }
                out << endl;
            }
        }
        else {
            // Legacy path: not reachable from the CLI (nextGen is always on) and
            // NOT compatible with the ordered parallel driver, since reduceAlleles
            // writes to cout directly. Kept only so the call site still compiles.
            var.reduceAlleles(
                varAlleles,
                variantFile,
                var,
                parseFlag,
                keepInfo,
                keepGeno,
                debug);
        }
        return true;
    }; // end process_record lambda

    // ---- Record-level parallel driver ----
    // Records are independent, but parsing mutates the shared file stream, so we
    // read a batch serially, decompose the expensive (realigned) records in the
    // batch concurrently (one thread per record; each WFA alignment runs
    // single-threaded), then emit results in input order. This scales ~linearly
    // with cores, unlike WFA-internal threading which plateaus around 2-3x. Peak
    // memory grows with the number of threads (each big alignment needs ~0.7GB).
    const int nthreads = max(1, thread_count);
    const size_t batch_size = (size_t)nthreads * 8;

    bool eof = false;
    while (!eof) {
        vector<WfaVariant> batch;
        batch.reserve(batch_size);
        for (size_t i = 0; i < batch_size; ++i) {
            WfaVariant v(variantFile);
            if (!variantFile.getNextVariant(v)) { eof = true; break; }
            batch.push_back(std::move(v));
        }
        if (batch.empty()) break;

        vector<string> outputs(batch.size());
        vector<string> errors(batch.size());
        vector<char> decomposed(batch.size(), 0); // 1 = realigned (output pre-rendered), 0 = pass-through

        #pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
        for (size_t i = 0; i < batch.size(); ++i) {
            // Exceptions must not escape an OpenMP region; capture and report
            // serially after the parallel loop.
            try {
                ostringstream oss;
                if (process_record(batch[i], oss)) {
                    decomposed[i] = 1;
                    outputs[i] = oss.str();
                }
            } catch (const std::exception& e) {
                errors[i] = e.what();
            }
        }

        // Emit in input order. Pass-through records are rendered here, straight to
        // cout, so cout's float precision carries across records exactly as in the
        // original serial code; a realigned record carries its own precision(2)
        // and then bumps cout to precision(2) for the records that follow it.
        for (size_t i = 0; i < batch.size(); ++i) {
            if (!errors[i].empty()) {
                cerr << "vcfwave: error decomposing record: " << errors[i] << endl;
                exit(1);
            }
            if (decomposed[i]) {
                cout << outputs[i];
                if (!outputs[i].empty()) cout.precision(2);
            } else {
                cout << batch[i] << endl;
            }
        }

        if (!quiet && file_size >= 0) {
            off_t pos = variantFile.file_pos();
            if (pos >= 0) {
                amount = (double)pos/(double)file_size;
                if (amount > prev_amount+0.003) {
                    prev_amount = amount;
                    print_progress(amount*100, start);
                }
            }
        }
    }

    if (!quiet) cerr << endl;

    return 0;
}
