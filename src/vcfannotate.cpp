#include "Variant.h"
#include "BedReader.h"
#include <getopt.h>

using namespace std;
using namespace vcflib;


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [<vcf file>]" << endl
         << endl
         << "options:" << endl 
         << "    -b, --bed   use annotations provided by this BED file" << endl
         << "    -k, --key   use this INFO field key for the annotations" << endl
         << "    -d, --default  use this INFO field key for records without annotations" << endl
         << endl
         << "Intersect the records in the VCF file with targets provided in a BED file." << endl
         << "Intersections are done on the reference sequences in the VCF file." << endl
         << "If no VCF filename is specified on the command line (last argument) the VCF" << endl
         << "read from stdin." << endl;
    exit(0);
}

int main(int argc, char** argv) {

    string bedFileName;
    string annotationInfoKey;
    string defaultAnnotationValue;

    if (argc == 1)
        printSummary(argv);

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"bed",  required_argument, 0, 'b'},
            {"key",  required_argument, 0, 'k'},
            {"default",  required_argument, 0, 'd'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hb:k:d:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
            case 'b':
                bedFileName = string(optarg);
                break;

            case 'k':
                annotationInfoKey = string(optarg);
                break;

            case 'd':
                defaultAnnotationValue = string(optarg);
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

    string line = "##INFO=<ID=" + annotationInfoKey + ",Number=1,Type=String,Description=\"Annotation from "
        + bedFileName + " delimited by ':'\">";
    variantFile.addHeaderLine(line);

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        BedTarget record(var.sequenceName, var.position, var.position + var.ref.size() - 1, "");
        vector<BedTarget*> overlaps = bed.targetsOverlapping(record);
        vector<string> annotations;
        if (!overlaps.empty()) {
            for (vector<BedTarget*>::iterator t = overlaps.begin(); t != overlaps.end(); ++t) {
                annotations.push_back((*t)->desc);
            }
            var.info[annotationInfoKey].push_back(join(annotations, ":"));
        } else if (!defaultAnnotationValue.empty()) {
            var.info[annotationInfoKey].push_back(defaultAnnotationValue);
        }
        cout << var << endl;
    }

    return 0;

}
