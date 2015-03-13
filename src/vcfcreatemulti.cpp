#include "Variant.h"
#include "convert.h"
#include <set>
#include <sstream>
#include <getopt.h>

using namespace std;
using namespace vcflib;


double convertStrDbl(const string& s) {
    double r;
    convert(s, r);
    return r;
}

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [file]" << endl
         << endl
         << "If overlapping alleles are represented across multiple records, merge" << endl
         << "them into a single record.  Currently only for indels." << endl;
    exit(0);
}

Variant createMultiallelic(vector<Variant>& vars) {

    if (vars.size() == 1) {
        return vars.front();
    }

    int maxpos = vars.front().position + vars.front().ref.size();
    for (vector<Variant>::iterator v = vars.begin(); v != vars.end(); ++v) {
        //cerr << *v << endl;
        if (maxpos < v->position + v->ref.size()) {
            maxpos = v->position + v->ref.size();
        }
    }

    int numalt = vars.size();
    //cerr << "gots overlapping vars " << vars.front().position << "-" << vars.back().position << endl;

    // get REF
    // use start position to extend all other alleles
    int start = vars.front().position;
    string ref = vars.front().ref;

    for (vector<Variant>::iterator v = vars.begin() + 1; v != vars.end(); ++v) {
        int sdiff = (v->position + v->ref.size()) - (start + ref.size());
        int pdiff = (start + ref.size()) - v->position;
        if (sdiff > 0) {
            ref.append(v->ref.substr(pdiff, sdiff));
        }
    }

    //cerr << "ref would be " << ref << " for vars from "
    //     << vars.front().position << " to " << vars.back().position << endl;

    Variant var = vars.front();
    var.alt.clear();
    var.ref = ref;

    for (vector<Variant>::iterator v = vars.begin(); v != vars.end(); ++v) {
        // add alternates and splice them into the reference
        int p5diff = v->position - var.position;
        int p3diff = (var.position + var.ref.size()) - (v->position + v->ref.size());
        string before;
        string after;
        if (p5diff > 0) {
            before = var.ref.substr(0, p5diff);
        }
        if (p3diff > 0 && p3diff < var.ref.size()) {
            after = var.ref.substr(var.ref.size() - p3diff);
        }
        if (p5diff || p3diff) {
            for (vector<string>::iterator a = v->alt.begin(); a != v->alt.end(); ++a) {
                var.alt.push_back(before);
                string& alt = var.alt.back();
                alt.append(*a);
                alt.append(after);
            }
        } else {
            for (vector<string>::iterator a = v->alt.begin(); a != v->alt.end(); ++a) {
                var.alt.push_back(*a);
            }
        }
    }

    stringstream s;
    s << vars.front().position << "-" << vars.back().position;
    var.info["combined"].push_back(s.str());

    return var;

}

int main(int argc, char** argv) {

    VariantCallFile variantFile;

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "h",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

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

    bool first = true;
    bool already = false;
    Variant var(variantFile);
    vector<Variant> vars;
    string lastSeq;

    while (variantFile.getNextVariant(var)) {

        if (lastSeq.empty()) {
            lastSeq = var.sequenceName;
        }

        if (vars.empty()) {
            vars.push_back(var);
            continue;
        } else {
            int maxpos = vars.front().position + vars.front().ref.size();
            for (vector<Variant>::iterator v = vars.begin(); v != vars.end(); ++v) {
                if (maxpos < v->position + v->ref.size()) {
                    maxpos = v->position + v->ref.size();
                }
            }
            if (var.sequenceName != lastSeq) {
                Variant result = createMultiallelic(vars);
                cout << result << endl;
                vars.clear();
                lastSeq = var.sequenceName;
                vars.push_back(var);
            } else if (var.position < maxpos) {
                vars.push_back(var);
            } else {
                Variant result = createMultiallelic(vars);
                cout << result << endl;
                vars.clear();
                vars.push_back(var);
            }
        }

    }

    if (!vars.empty()) {
        Variant result = createMultiallelic(vars);
        cout << result << endl;
    }

    return 0;

}

