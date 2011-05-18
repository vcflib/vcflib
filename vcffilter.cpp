#include "Variant.h"
#include "split.h"
#include <getopt.h>
#include <deque>

using namespace std;
using namespace vcf;

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "options:" << endl 
         << "    -f, --info-filter     specifies a filter to apply to the info fields of records" << endl
         << "    -g, --genotype-filter specifies a filter to apply to the genotype fields of records" << endl
         << "    -t, --tag             tag vcf records as filtered with this tag instead of suppressing them" << endl
         << "    -v, --invert          inverts the filter, e.g. grep -v" << endl
         << "    -o, --or              use logical OR instead of AND to combine filters" << endl
         << endl
         << "Filter the specified vcf file using the set of filters." << endl
         << "Filters are specified in the form \"<ID> <operator> <value>:" << endl
         << " -f \"DP > 10\"  # for info fields" << endl
         << " -g \"GT = 1|1\" # for genotype fields" << endl
         << " -f \"CpG\"  # for 'flag' fields" << endl
         << endl
         << "Operators can be any of: =, !, <, >, |, &" << endl
         << endl
         << "Any number of filters may be specified.  They are combined via logical AND" << endl
         << "unless --or is specified on the command line.  Obtain logical negation through" << endl
         << "the use of parentheses, e.g. !( DP = 10)." << endl
         << endl
         << "For convenience, you can specify \"QUAL\" to refer to the quality of the site, even" << endl
         << "though it does not appear in the INFO fields." << endl
         << endl;
    exit(0);
}

bool passesFilters(Variant& var, vector<VariantFilter>& filters, bool logicalOr) {
    for (vector<VariantFilter>::iterator f = filters.begin(); f != filters.end(); ++f) {
        string s = "";
        if (logicalOr) {
            if (f->passes(var, s)) return true;
        } else {
            if (!f->passes(var, s)) return false;
        }
    }
    if (logicalOr)
        return false;
    else
        return true;
}


int main(int argc, char** argv) {

    int c;
    bool invert = false;
    bool logicalOr = false;
    vector<string> infofilterStrs;
    vector<VariantFilter> infofilters;
    vector<string> genofilterStrs;
    vector<VariantFilter> genofilters;
    string tag = "";
    string filterSpec;

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
            {"or", no_argument, 0, 'o'},
            //{"length",  no_argument, &printLength, true},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvof:g:t:",
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
            filterSpec += " " + string(optarg);
            infofilterStrs.push_back(string(optarg));
            break;
 
          case 'g':
            filterSpec += " genotypes filtered with: " + string(optarg);
            genofilterStrs.push_back(string(optarg));
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

          case 'o':
            logicalOr = true;
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

    filterSpec = filterSpec.substr(1); // strip leading " "

    VariantCallFile variantFile;
    if (optind == argc - 1) {
        string inputFilename = argv[optind];
        variantFile.open(inputFilename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    for (vector<string>::iterator f = infofilterStrs.begin(); f != infofilterStrs.end(); ++f) {
        infofilters.push_back(VariantFilter(*f, VariantFilter::RECORD, variantFile.infoTypes));
    }

    for (vector<string>::iterator f = genofilterStrs.begin(); f != genofilterStrs.end(); ++f) {
        genofilters.push_back(VariantFilter(*f, VariantFilter::SAMPLE, variantFile.formatTypes));
    }

    vector<string> headerlines = split(variantFile.header, "\n");
    variantFile.header.clear();
    for (vector<string>::iterator l = headerlines.begin(); l != headerlines.end(); ++l) {
        if (!filterSpec.empty() && (l->find("INFO") != string::npos || l + 1 == headerlines.end())) {
            variantFile.header += "##filter=\"" + filterSpec + "\"\n";
            filterSpec.clear();
        }
        variantFile.header += *l + ((l + 1 == headerlines.end()) ? "" : "\n");
    }
    cout << variantFile.header << endl;

    // holds records from the same site--- output the records chunk-wise
    deque<string> variantOutputDeque;

    int lastPosition = 0;
    string lastSequenceName;
    bool lastPositionPassed = true;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        //cout << var.position << " " << lastPosition << endl;
        if (var.position > lastPosition || var.sequenceName != lastSequenceName) {
            if (lastPositionPassed) {
                while (!variantOutputDeque.empty()) {
                    string& v = variantOutputDeque.back();
                    cout << v << endl;
                    variantOutputDeque.pop_back();
                }
            } else {
                variantOutputDeque.clear();
            }
            lastPositionPassed = true;
            lastPosition = var.position;
            lastSequenceName = var.sequenceName;
        }
        if (!genofilters.empty()) {
            for (vector<VariantFilter>::iterator f = genofilters.begin(); f != genofilters.end(); ++f) {
                f->removeFilteredGenotypes(var);
            }
        }
        if (!infofilters.empty()) {
            bool passes = passesFilters(var, infofilters, logicalOr);
            if (invert) {
                passes = !passes;
            }
            if (passes) {
                if (!tag.empty()) {
                    var.addFilter(tag);
                    stringstream s; s << var;
                    variantOutputDeque.push_front(s.str());
                } else {
                    variantOutputDeque.push_front(variantFile.line);
                }
            } else if (!tag.empty()) {
                variantOutputDeque.push_front(variantFile.line);
            }
            lastPositionPassed &= passes;
        } else {
            if (genofilters.empty()) {
                variantOutputDeque.push_front(variantFile.line);
            } else {
                stringstream s; s << var;
                variantOutputDeque.push_front(s.str());
            }
        }
    }

    if (lastPositionPassed) {
        while (!variantOutputDeque.empty()) {
            string& v = variantOutputDeque.back();
            cout << v << endl;
            variantOutputDeque.pop_back();
        }
    }

    return 0;

}

