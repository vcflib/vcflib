/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020-2023 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
extern "C" {
#include "vcf-c-api.h"
}
#include <set>
#include <sstream>
#include <getopt.h>
#include "progress.h"

using namespace std;
using namespace vcflib;

#ifdef NO_ZIG
bool nextGen  = false;
#else
bool nextGen  = true;
#endif

bool quiet = false;
off_t file_size = -1;

void printSummary(char** argv) {
    cerr << R"(
Usage: vcfcreatemulti [options] [file]

Go through sorted VCF and when overlapping alleles are represented across multiple records, merge them into a single multi-ALT record. See the documentation for more information.

options:

    --quiet           no progress bar
    --legacy          legacy mode (old C++ implementation does not do genotypes)

Type: transformation
)";
    exit(1);
}

// Helper function to convert a vector of objects to a vector of
// pointers (for zig).
vector<void *> ptr_vec(vector<Variant> &vars) {
    vector<void *> ptrs;
    for (auto &v: vars) {
        ptrs.push_back(&v);
    }
    return ptrs;
}

Variant createMultiallelic_legacy(vector<Variant>& vars) {

    if (vars.size() == 1) {
        return vars.front();
    }

    Variant first = vars.front();
    Variant nvar = first;

    // get REF
    // use start position to extend all other alleles
    int start = first.position;
    string ref = first.ref;

    // expand reference using all references in list
    for (vector<Variant>::iterator v = vars.begin() + 1; v != vars.end(); ++v) {
        int sdiff = (v->position + v->ref.size()) - (start + ref.size());
        int pdiff = (start + ref.size()) - v->position;
        if (sdiff > 0) {
            ref.append(v->ref.substr(pdiff, sdiff));
        }
    }

    Variant mvar = first;
    mvar.alt.clear();
    mvar.ref = ref;

    // Correct alts using the new reference
    for (const auto& v : vars) {
        // add alternates and splice them into the reference
        int p5diff = v.position - mvar.position;
        int p3diff = (mvar.position + mvar.ref.size()) - (v.position + v.ref.size());
        string before;
        string after;
        if (p5diff > 0) {
            before = mvar.ref.substr(0, p5diff);
        }
        if (p3diff > 0 && p3diff < mvar.ref.size()) {
            after = mvar.ref.substr(mvar.ref.size() - p3diff);
        }
        if (p5diff || p3diff) {
            for (const auto& a : v.alt) {
                mvar.alt.push_back(before);
                string& alt = mvar.alt.back();
                alt.append(a);
                alt.append(after);
            }
        } else {
            for (const auto& a : v.alt) {
                mvar.alt.push_back(a);
            }
        }
    }

    stringstream s;
    s << vars.front().position << "-" << vars.back().position;
    mvar.info["combined"].push_back(s.str());

    return mvar;
}

#ifndef NO_ZIG
Variant createMultiallelic_zig(vector<Variant>& vars) {

    if (vars.size() == 1) {
        return vars.front();
    }
    Variant first = vars.front();
    Variant nvar = first;

    Variant *mvar = (Variant *)zig_create_multi_allelic(&nvar, ptr_vec(vars).data() , vars.size());
    stringstream s;
    s << vars.front().position << "-" << vars.back().position;
    mvar->info["combined"].push_back(s.str());

    return *mvar;
}
#endif

// This function is called for every line/variant in the VCF file
Variant createMultiallelic(vector<Variant>& vars) {
    #ifndef NO_ZIG
    if (nextGen)
        return createMultiallelic_zig(vars);
    else
    #endif
        return createMultiallelic_legacy(vars);
}

int main(int argc, char** argv) {

    VariantCallFile variantFile;

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"quiet", no_argument, 0, 'q'},
            {"legacy", no_argument, 0, 'l'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hnq",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

            case 'l':
                nextGen = false;
                break;

            case 'q':
                quiet = true;
                break;

            case 'h':
                printSummary(argv);
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
        file_size = get_file_size(filename.c_str());
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    variantFile.addHeaderLine("##INFO=<ID=combined,Number=1,Type=String,Description=\"Range of overlapping variants which were combined into this one using vcfcreatemulti.\">");
    variantFile.addHeaderLine("##INFO=<ID=MULTI,Number=1,Type=String,Description=\"Identify problematic ALT reconstruction - see vcfcreatemulti.md doc.\">");

    cout << variantFile.header << endl;

    Variant var(variantFile);
    string lastSeqName;
    vector<Variant> vars;

    double amount = 0.0, prev_amount = 0.0;
    uint64_t start = get_timestamp();
    string prev_chr = "none";
    size_t prev_pos = 0;

    if (!quiet)
        cerr << "vcfcreatemulti " << VCFLIB_VERSION << " processing..." << endl;

    while (variantFile.getNextVariant(var)) {

        if (prev_pos && prev_chr == var.sequenceName && prev_pos > var.position) {
            cerr << "ERROR: VCF data is not sorted! at " << prev_chr << ":" << prev_pos << " and " << var.sequenceName << ":" << var.position << endl;
            exit(8);
        }

        prev_chr = var.sequenceName;
        prev_pos = var.position;

        amount = (double)variantFile.file_pos()/(double)file_size;
        // cerr << file_size << "," << variantFile.file_pos() << "=" << amount << endl;
        if (!quiet && variantFile.file_pos() >= 0 && file_size >= 0 && amount > prev_amount+0.003) {
            prev_amount = amount;
            print_progress(amount*100, start);
        }

        if (lastSeqName.empty()) { // track the previous sequence name
            lastSeqName = var.sequenceName;
        }

        if (vars.empty()) { // track list of variants in window (alt alleles)
            vars.push_back(var);
            Variant *copy2 = &vars.back();
            continue;
        } else {
            // compute maxpos as the most right position in the current reference window. Note
            // that the window may get expanded at every step.
            auto first = vars.front();
            auto maxpos = first.position + first.ref.size();
            for (const auto& v: vars) {
                if (maxpos < v.position + v.ref.size()) {
                    maxpos = v.position + v.ref.size();
                }
            }
            if (var.sequenceName != lastSeqName) {
                // next chromosome, contig or whatever
                Variant result = createMultiallelic(vars);
                cout << result << endl;
                vars.clear();
                lastSeqName = var.sequenceName;
                vars.push_back(var);
                Variant *copy2 = &vars.back();
            } else if (var.position < maxpos) {
                // As long as it is in window add it to the list
                vars.push_back(var);
                Variant *copy2 = &vars.back();
            } else {
                // Next variant is out of window, so create single line variant
                Variant result = createMultiallelic(vars);
                cout << result << endl;
                vars.clear();
                vars.push_back(var);
                Variant *copy2 = &vars.back();
            }
        }
    }

    if (!vars.empty()) {
        Variant result = createMultiallelic(vars);
        cout << result << endl;
    }

    #ifndef NO_ZIG
    if (nextGen) {
        zig_display_warnings();
        zig_cleanup();
    }
    #endif

    if (!quiet) cerr << endl;

    return 0;

}
