#include "Variant.h"
#include "BedReader.h"
#include <getopt.h>
#include "mt19937ar.h"
#include <sstream>
#include <iostream>
#include "convert.h"

using namespace std;
using namespace vcflib;


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [<vcf file>]" << endl
         << endl
         << "options:" << endl 
         << "    -r, --rate RATE          base sampling probability per locus" << endl
         << "    -s, --scale-by KEY       scale sampling likelihood by this Float info field" << endl
         << "    -p, --random-seed N      use this random seed (by default read from /dev/random)" << endl
         << "    -q, --pseudorandom-seed  use a pseudorandom seed (by default read from /dev/random)" << endl
         << endl
         << "Randomly sample sites from an input VCF file, which may be provided as stdin." << endl
         << "Scale the sampling probability by the field specified in KEY.  This may be" << endl
         << "used to provide uniform sampling across allele frequencies, for instance." << endl;
    exit(0);
}

int main(int argc, char** argv) {

    double rate = 1.0;
    int seed = 0;
    bool useprng = false;
    string scaleByKey;

    if (argc == 1)
        printSummary(argv);

    int c;
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"rate",  required_argument, 0, 'r'},
                {"scale-by",  required_argument, 0, 's'},
                {"random-seed",  required_argument, 0, 'p'},
                {"pseudorandom-seed",  required_argument, 0, 'q'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hqr:s:p:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
        case 'r':
            rate = atof(optarg);
            break;

        case 's':
            scaleByKey = optarg;
            break;

        case 'p':
            seed = atoi(optarg);
            break;

        case 'q':
            useprng = true;
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

    // seed prng with random bits from /dev/random
    if (!seed) {
        fstream random;
        if (useprng) {
            random.open("/dev/urandom", fstream::in);
        } else {
            random.open("/dev/random", fstream::in);
        }
        random.get((char*) &seed, sizeof(int));
        random.close();
    }

    init_genrand(seed);

    vector<string> args;
    for (int i = 0; i < argc; ++i) {
        args.push_back(argv[i]);
    }

    stringstream liness;
    liness << "##sampling=\"random sampling using "
           << join(args, " ")
           << " using random seed "
           << seed << "\"";
    variantFile.addHeaderLine(liness.str());

    cout << variantFile.header << endl;
    
    // check that we can use the scaling key
    if (!scaleByKey.empty()) {
        if (variantFile.infoTypes.find(scaleByKey) == variantFile.infoTypes.end()) {
            cerr << "could not find info key " << scaleByKey << endl;
            exit(1);
        } else {
            if (variantFile.infoTypes[scaleByKey] != FIELD_FLOAT) {
                cerr << "cannot use " << scaleByKey << " as a scaling factor, as it is not of type Float" << endl;
                exit(1);
            }
        }
    }

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        double randN = genrand_real1();
        if (!scaleByKey.empty()) {
            if (var.info.find(scaleByKey) != var.info.end()) {
                double val;

                // hack, sum the values of interest if we have multiple values
                // really, this is only suitable for AF stuff
                vector<string>& vals = var.info[scaleByKey];
                for (vector<string>::iterator b = vals.begin(); b != vals.end(); ++b) {
                    double f;
                    convert(*b, f);
                    val += f;
                }
                val /= vals.size();

                if (val > 1) {
                    cerr << "cannot scale by " << scaleByKey << "=" << val << " as it is > 1" << endl;
                    exit(1);
                }
                randN *= val;
            }
        }
        if (randN < rate) {
            cout << var << endl;
        }
    }

    return 0;

}
