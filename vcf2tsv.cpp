#include "Variant.h"
#include <getopt.h>

using namespace std;
using namespace vcf;

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [-n null_string]" << " [vcf file]" << endl
	 << "Converts stdin or given VCF file to tab-delimited format, using null string to replace empty values in the table." << endl;
    exit(1);
}


int main(int argc, char** argv) {

    string nullval;

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
	    {"null-value", required_argument, 0, 'n'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hn:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

	    case 'n':
	        nullval = optarg;
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

    VariantCallFile variantFile;
    bool usingstdin = false;
    string inputFilename;
    if (optind == argc - 1) {
        inputFilename = argv[optind];
        variantFile.open(inputFilename);
    } else {
        if (!variantFile.open(std::cin)) {
	    if (argc == 1) {
		printSummary(argv);
	    } else {
		cerr << "could not open stdin for reading as VCF" << endl;
		exit(1);
	    }
	}
	usingstdin = true;
    }

    if (!variantFile.is_open()) {
        return 1;
    }
    // obtain all possible field names
    vector<string> infofields;
    vector<string> infoflags;

    for (map<string, VariantFieldType>::iterator i = variantFile.infoTypes.begin(); i != variantFile.infoTypes.end(); ++i) {
        if (i->second == FIELD_BOOL) {
            infoflags.push_back(i->first);
        } else {
            infofields.push_back(i->first);
        }
    }

    // write header

    // defaults
    cout << "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER";
    
    // configurable info field
    for (vector<string>::iterator i = infofields.begin(); i != infofields.end(); ++i) {
        cout << "\t" << *i;
    }
    for (vector<string>::iterator i = infoflags.begin(); i != infoflags.end(); ++i) {
        cout << "\t" << *i;
    }
    cout << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

        int altindex = 0;
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a, ++altindex) {

            string& altallele = *a;

            cout << var.sequenceName << "\t"
                 << var.position << "\t"
                 << var.id << "\t"
                 << var.ref << "\t"
                 << altallele << "\t"
                 << var.quality << "\t"
                 << var.filter;

            for (vector<string>::iterator i = infofields.begin(); i != infofields.end(); ++i) {
                vector<string> value;
                string& name = *i;
                map<string, vector<string> >::iterator f = var.info.find(name);
                if (f != var.info.end()) {
                    value = f->second;
                    if (value.size() == 1) {
                        cout << "\t" << value.front();
                    } else if (value.size() == var.alt.size()) {
                        cout << "\t" << value.at(altindex);
                    } else {
                        cout << "\t" << nullval; // null
                    }
                } else {
                    cout << "\t" << nullval; // null
                }
            }

            for (vector<string>::iterator i = infoflags.begin(); i != infoflags.end(); ++i) {
                string value;
                string& name = *i;
                map<string, bool>::iterator f = var.infoFlags.find(name);
                cout << "\t";
                if (f != var.infoFlags.end()) {
                    cout << 1;
                } else {
                    cout << 0;
                }
            }

            cout << endl;

        }

    }

    return 0;

}

