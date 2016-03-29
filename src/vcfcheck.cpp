#include "Variant.h"
#include "split.h"
#include "Fasta.h"
#include "gpatInfo.hpp"
#include <getopt.h>

using namespace std;
using namespace vcflib;

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "options:" << endl 
         << "    -f, --fasta-reference  FASTA reference file to use to obtain primer sequences" << endl
         << "    -x, --exclude-failures If a record fails, don't print it.  Otherwise do." << endl
         << "    -k, --keep-failures    Print if the record fails, otherwise not." << endl
	 << "    -h, --help       Print this message." << endl
	 << "    -v, --version    Print version." << endl
         << endl
         << "Verifies that the VCF REF field matches the reference as described." << endl
         << endl;
    exit(0);
}


int main(int argc, char** argv) {

    int c;
    string fastaRef;
    bool keepFailures = false;
    bool excludeFailures = false;

    if (argc == 1)
        printSummary(argv);

    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"fasta-reference",  required_argument, 0, 'f'},
                {"exclude-failures",  no_argument, 0, 'x'},
                {"keep-failures",  no_argument, 0, 'k'},
                {"version",  no_argument, 0, 'v'},
                //{"length",  no_argument, &printLength, true},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvxkf:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
            printf ("\n");
            break;

        case 'f':
	  {
            fastaRef = optarg;
            break;
	  }
	case 'v':
	  {
	    printBasicVersion();
	    exit(0);
	  }
        case 'x':
	  {
            excludeFailures = true;
            break;
	  }
        case 'k':
	  {
            keepFailures = true;
            break;
	  }
        case 'h':
	  {
            printSummary(argv);
            exit(0);
            break;
	  }
        case '?':
	  {
            /* getopt_long already printed an error message. */
            printSummary(argv);
            exit(1);
            break;
	  }
        default:
            abort ();
        }
    }

    if (fastaRef.empty()) {
        cerr << "a FASTA reference sequence must be specified" << endl;
        exit(1);
    }

    FastaReference ref;
    ref.open(fastaRef);

    VariantCallFile variantFile;
    string inputFilename;
    if (optind == argc - 1) {
        inputFilename = argv[optind];
        variantFile.open(inputFilename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    if (keepFailures || excludeFailures) {
        cout << variantFile.header << endl;
    }

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        int refstart = var.position - 1; // convert to 0-based
        string matchedRef = ref.getSubSequence(var.sequenceName, refstart, var.ref.size());
        if (var.ref != matchedRef) {
            if (keepFailures) {
                cout << var << endl;
            } else if (!excludeFailures) {
                cout << "mismatched reference " << var.ref << " should be " << matchedRef << " at "
                     << var.sequenceName << ":" << var.position << endl;
            }
        } else if (excludeFailures) {
            cout << var << endl;
        }
    }

    return 0;

}

