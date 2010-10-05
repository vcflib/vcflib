#include "Variant.h"
#include "Split.h"
#include <getopt.h>


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "options:" << endl 
         << "    -f, --info-filter     specifies a filter to apply to the info fields of records" << endl
         << "    -g, --genotype-filter specifies a filter to apply to the genotype fields of records" << endl
         << "    -t, --tag             tag vcf records as filtered with this tag instead of suppressing them" << endl
         << "    -v, --invert          inverts the filter, e.g. grep -v" << endl
         << endl
         << "Filter the specified vcf file using the set of filters." << endl
         << endl;
    exit(0);
}

bool passesFilters(Variant& var, vector<VariantFilter>& filters) {
    for (vector<VariantFilter>::iterator f = filters.begin(); f != filters.end(); ++f) {
        if (!f->passes(var))
            return false;
    }
    return true;
}


int main(int argc, char** argv) {

    int c;
    bool invert = false;
    vector<VariantFilter> infofilters;
    vector<VariantFilter> genofilters;
    string tag = "";

    if (argc == 1)
        printSummary(argv);

    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"info-filter",  required_argument, 0, 'f'},
            {"genotype-filter",  required_argument, 0, 'g'},
            {"tag", required_argument, 0, 't'},
            {"invert", no_argument, 0, 'v'},
            //{"length",  no_argument, &printLength, true},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvf:g:t:",
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
            infofilters.push_back(VariantFilter(optarg, true));
            break;
 
          case 'g':
            genofilters.push_back(VariantFilter(optarg, false));
            break;
 
          case 't':
            tag = optarg;
            break;
 
          case 'h':
            printSummary(argv);
            exit(0);
            break;

          case 'v':
            invert = true;
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

    string inputFilename;
    if (optind == argc - 1) {
        inputFilename = argv[optind];
    } else {
        cerr << "No input file provided." << endl;
        exit(1);
    }

    VariantCallFile variantFile;
    if (!variantFile.openVCF(inputFilename)) {
        return 1;
    }

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        bool passesInfo = passesFilters(var, infofilters);
        bool passesGeno = passesFilters(var, genofilters);
        bool passes = (passesInfo && passesGeno);
        if (invert) {
            passes = !passes;
        }
        if (passes) {
            var.addFilter(tag);
            cout << var << endl;
        }
    }

    return 0;

}

