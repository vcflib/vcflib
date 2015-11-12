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
         << "    -f, --fasta-reference  FASTA reference file to use to obtain primer sequences" << endl
         << "    -l, --primer-length    The length of the primer sequences on each side of the variant" << endl
         << endl
         << "For each VCF record, extract the flanking sequences, and write them to stdout as FASTA" << endl
         << "records suitable for alignment.  This tool is intended for use in designing validation" << endl
         << "experiments.  Primers extracted which would flank all of the alleles at multi-allelic" << endl
         << "sites.  The name of the FASTA \"reads\" indicates the VCF record which they apply to." << endl
         << "The form is >CHROM_POS_LEFT for the 3' primer and >CHROM_POS_RIGHT for the 5' primer," << endl
         << "for example:" << endl
         << endl
         << ">20_233255_LEFT" << endl
         << "CCATTGTATATATAGACCATAATTTCTTTATCCAATCATCTGTTGATGGA" << endl
         << ">20_233255_RIGHT" << endl
         << "ACTCAGTTGATTCCATACCTTTGCCATCATGAATCATGTTGTAATAAACA" << endl
         << endl;
    exit(0);
}


int main(int argc, char** argv) {

    int c;
    string fastaRef;
    int primerLength = 0;

    if (argc == 1)
        printSummary(argv);

    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"fasta-reference",  required_argument, 0, 'f'},
            {"primer-length", required_argument, 0, 'l'},
            //{"length",  no_argument, &printLength, true},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hf:l:",
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

          case 'l':
            primerLength = atoi(optarg);
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

    if (primerLength == 0) {
        cerr << "a primer length must be specified" << endl;
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

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        // get the ref start and end positions
        int refstart = var.position - 1; // convert to 0-based
        int refend = var.position + var.ref.size() - 1;
        string leftprimer = ref.getSubSequence(var.sequenceName, refstart - primerLength, primerLength);
        string rightprimer = ref.getSubSequence(var.sequenceName, refend, primerLength);
        //cout << var << endl;
        cout << ">" << var.sequenceName << "_" << var.position << "_LEFT" << endl
             << leftprimer << endl
             << ">" << var.sequenceName << "_" << var.position << "_RIGHT" << endl
             << rightprimer << endl;
    }

    return 0;

}

