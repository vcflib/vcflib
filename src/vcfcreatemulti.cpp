/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020-2022 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
extern "C" {
#include "vcf-c-api.h"
}
#include <set>
#include <sstream>
#include <getopt.h>

using namespace std;
using namespace vcflib;

bool nextGen  = false;

// extern "C" void *zig_create_multi_allelic(Variant *retvar, Variant *varlist[], long size);

void printSummary(char** argv) {
    cerr << R"(
Usage: vcfcreatemulti [options] [file]

Go through sorted VCF and if overlapping alleles are represented
across multiple records, merge them into a single record.  Currently
only for indels.

options:
     -n, --nextgen           next gen mode.

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
    for (vector<Variant>::iterator v = vars.begin(); v != vars.end(); ++v) {
        // add alternates and splice them into the reference
        int p5diff = v->position - mvar.position;
        int p3diff = (mvar.position + mvar.ref.size()) - (v->position + v->ref.size());
        string before;
        string after;
        if (p5diff > 0) {
            before = mvar.ref.substr(0, p5diff);
        }
        if (p3diff > 0 && p3diff < mvar.ref.size()) {
            after = mvar.ref.substr(mvar.ref.size() - p3diff);
        }
        if (p5diff || p3diff) {
            for (vector<string>::iterator a = v->alt.begin(); a != v->alt.end(); ++a) {
                mvar.alt.push_back(before);
                string& alt = mvar.alt.back();
                alt.append(*a);
                alt.append(after);
            }
        } else {
            for (vector<string>::iterator a = v->alt.begin(); a != v->alt.end(); ++a) {
                mvar.alt.push_back(*a);
            }
        }
    }

    stringstream s;
    s << vars.front().position << "-" << vars.back().position;
    mvar.info["combined"].push_back(s.str());

    return mvar;
}

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


Variant createMultiallelic(vector<Variant>& vars) {
    if (nextGen)
        return createMultiallelic_zig(vars);
    else
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
            {"nextgen", no_argument, 0, 'n'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hn",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

            case 'n':
                nextGen = true;
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
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    variantFile.addHeaderLine("##INFO=<ID=combined,Number=1,Type=String,Description=\"Range of overlapping variants which were combined into this one using vcfcreatemulti.\">");

    cout << variantFile.header << endl;

    Variant var(variantFile);
    string lastSeqName;
    vector<Variant> vars;

    // cout << "Calling into Zig" << endl;
    //auto var_window = zig_variant_window();

    while (variantFile.getNextVariant(var)) {

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
            for (auto v: vars) {
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

    return 0;

}
