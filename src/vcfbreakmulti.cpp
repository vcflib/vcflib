#include "Variant.h"
#include "convert.h"
#include <set>
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
	 << "If multiple alleles are specified in a single record, break the record into" << endl
	 << "multiple lines, preserving allele-specific INFO fields." << endl;
    exit(0);
}

int main(int argc, char** argv) {

    bool includePreviousBaseForIndels = true;
    bool useMNPs = false;

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

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

        int numalt = var.alt.size();

        if (numalt == 1) {
            cout << var << endl;
            continue;
        }

        vector<Variant> variants;
        for (int i = 0; i < numalt; ++i) {
            variants.push_back(var);
        }

        for (int i = 0; i < numalt; ++i) {
            Variant& v = variants.at(i);
            vector<string> altsToRemove;
            for (int j = 0; j < numalt; ++j) {
                if (j != i) {
                    altsToRemove.push_back(var.alt.at(j));
                }
            }
            for (vector<string>::iterator a = altsToRemove.begin(); a != altsToRemove.end(); ++a) {
                v.removeAlt(*a);
            }
        }

        for (vector<Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
            cout << *v << endl;
        }
    }

    return 0;

}

