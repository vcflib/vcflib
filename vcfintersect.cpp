#include "Variant.h"
#include "BedReader.h"
#include <getopt.h>

using namespace std;
using namespace vcf;


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [<vcf file>]" << endl
         << endl
         << "options:" << endl 
         << "    -b, --bed          use intervals provided by this BED file" << endl
         << "    -i, --inverse      invert the selection, printing only records which would" << endl
         << "                       not have been printed out" << endl
         << endl
         << "Intersect the records in the VCF file with targets provided in a BED file." << endl
         << "Intersections are done on the reference sequences in the VCF file." << endl
         << "If no VCF filename is specified on the command line (last argument) the VCF" << endl
         << "read from stdin." << endl;
    exit(0);
}

int main(int argc, char** argv) {

    string bedFileName;
    bool invert = false;
    bool contained = true;
    bool overlapping = false;

    if (argc == 1)
        printSummary(argv);

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"invert", no_argument, 0, 'i'},
            {"bed",  required_argument, 0, 'b'},
            {"vcf",  required_argument, 0, 'v'},
            {"contained",  no_argument, 0, 'c'},
            {"overlapping", no_argument, 0, 'o'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvcob:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
            case 'b':
                bedFileName = string(optarg);
                break;

            case 'v':
                invert = true;
                break;

            case 'c':
                contained = true;
                break;

            case 'o':
                overlapping = true;
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

    if (bedFileName.empty()) {
        cerr << "a BED file is required when intersecting" << endl;
        exit(1);
    }

    BedReader bed(bedFileName);

    VariantCallFile variantFile;
    string inputFilename;
    if (optind == argc - 1) {
        inputFilename = argv[optind];
        variantFile.open(inputFilename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        cout << "could not open VCF file" << endl;
        return 1;
    }

    variantFile.parseSamples = false; // faster

    cout << variantFile.header;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        BedTarget record(var.sequenceName, var.position, var.position + var.ref.size(), "");
        vector<BedTarget*> overlaps = bed.targetsOverlapping(record);
        if (!invert && !overlaps.empty()) {
            cout << variantFile.line << endl;
        } else if (invert && overlaps.empty()) {
            cout << variantFile.line << endl;
        }
    }

    return 0;

}

