#include "Variant.h"
#include "split.h"
#include "Fasta.h"
#include <getopt.h>

using namespace std;
using namespace vcflib;

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "options:" << endl 
         << "    -f, --fasta-reference  FASTA reference file to use to obtain flanking sequences" << endl
         << "    -w, --window-size      Size of the window over which to calculate entropy" << endl
         << endl
         << "Anotates the output VCF file with, for each record, EntropyLeft, EntropyRight," << endl
         << "EntropyCenter, which are the entropies of the sequence of the given window size to the" << endl
         << "left, right, and center  of the record.  Also adds EntropyRef and EntropyAlt for each alt." << endl
         << endl;
    exit(0);
}


int main(int argc, char** argv) {

    int c;
    string fastaRef;
    int windowSize = 0;

    if (argc == 1)
        printSummary(argv);

    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"fasta-reference",  required_argument, 0, 'f'},
            {"window-size", required_argument, 0, 'w'},
            //{"length",  no_argument, &printLength, true},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hf:w:",
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
            fastaRef = optarg;
            break;

          case 'w':
            windowSize = atoi(optarg);
            break;
 
          case 'h':
            printSummary(argv);
            exit(0);
            break;

          case '?':
            /* getopt_long already printed an error message. */
            printSummary(argv);
            exit(1);
            break;
 
          default:
            abort ();
          }
      }

    if (windowSize == 0) {
        cerr << "a window size must be specified" << endl;
        exit(1);
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

    variantFile.addHeaderLine("##INFO=<ID=EntropyLeft,Number=1,Type=Float,Description=\"Entropy of left-flanking sequence of "+ convert(windowSize) +"bp\">");
    variantFile.addHeaderLine("##INFO=<ID=EntropyCenter,Number=1,Type=Float,Description=\"Entropy of centered sequence of "+ convert(windowSize) +"bp\">");
    variantFile.addHeaderLine("##INFO=<ID=EntropyRight,Number=1,Type=Float,Description=\"Entropy of right-flanking sequence of "+ convert(windowSize) +"bp\">");
    variantFile.addHeaderLine("##INFO=<ID=EntropyRef,Number=1,Type=Float,Description=\"Entropy of REF allele\">");
    variantFile.addHeaderLine("##INFO=<ID=EntropyAlt,Number=A,Type=Float,Description=\"Entropy of each ALT allele\">");

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

        // get the ref start and end positions
        int refstart = var.position - 1; // convert to 0-based
        int refend = var.position + var.ref.size() - 1;
        string leftseq = ref.getSubSequence(var.sequenceName, refstart - windowSize, windowSize);
        string rightseq = ref.getSubSequence(var.sequenceName, refend, windowSize);
        string centerseq = ref.getSubSequence(var.sequenceName, refstart - windowSize/2, windowSize);
        double entropyLeft = shannon_H((char*) &leftseq[0], windowSize);
        double entropyRight = shannon_H((char*) &rightseq[0], windowSize);
        double entropyCenter = shannon_H((char*) &centerseq[0], windowSize);
        double entropyRef = shannon_H((char*) var.ref.c_str(), var.ref.size());

        var.info["EntropyLeft"].clear();
        var.info["EntropyRight"].clear();
        var.info["EntropyCenter"].clear();
        var.info["EntropyRef"].clear();
        var.info["EntropyAlt"].clear();

        var.info["EntropyLeft"].push_back(convert(entropyLeft));
        var.info["EntropyRight"].push_back(convert(entropyRight));
        var.info["EntropyCenter"].push_back(convert(entropyCenter));
        var.info["EntropyRef"].push_back(convert(entropyRef));

        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            double entropyAlt = shannon_H((char*) a->c_str(), a->size());
            var.info["EntropyAlt"].push_back(convert(entropyAlt));
        }

        cout << var << endl;
    }

    return 0;

}

