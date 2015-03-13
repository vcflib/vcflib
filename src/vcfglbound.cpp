#include "Variant.h"
#include "split.h"
#include <getopt.h>

using namespace std;
using namespace vcflib;

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "options:" << endl 
         << "    -b, --bound N          Bound GLs to this limit." << endl
         << "    -x, --exclude-broken   If GLs are > 0, remove site." << endl
         << endl
         << "Adjust GLs so that the maximum GL is 0 by dividing all GLs for each sample by the max." << endl
         << "Then cap (bound) at N (e.g. -10)." << endl;
    exit(0);
}


int main(int argc, char** argv) {

    bool excludeBroken = false;
    double glBound = 0;
    int c;

    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"bound",  required_argument, 0, 'b'},
            {"exclude-broken",  no_argument, 0, 'x'},
            //{"length",  no_argument, &printLength, true},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hxb:",
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

          case 'b':
              glBound = atof(optarg);
              break;
          
          case 'x':
              excludeBroken = true;
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

    if (glBound == 0) {
        cerr << "a bound is required when running vcfglbound (try -10)" << endl;
        exit(1);
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
        return 1;
    }

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        if (find(var.format.begin(), var.format.end(), "GL") == var.format.end()) {
            cout << var << endl;
            continue;
        }
        if (find(var.format.begin(), var.format.end(), "GT") == var.format.end()) {
            var.format.push_back("GT");
            reverse(var.format.begin(), var.format.end());
        }
        bool isbroken = false;
        for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin();
             s != var.samples.end(); ++s) {
            map<string, vector<string> >& sample = s->second;
            map<string, vector<string> >::iterator l = sample.find("GL");
            if (l != sample.end()) {

                // find the gl max
                vector<string>& glstrs = l->second;
                vector<double> gls;
                for (vector<string>::iterator gl = glstrs.begin(); gl != glstrs.end(); ++gl) {
                    double d;
                    convert(*gl, d);
                    gls.push_back(d);
                }

                isbroken = false; // reset every iteration
                for (vector<double>::iterator g = gls.begin(); g != gls.end(); ++g) {
                    if (*g > 0) {
                        isbroken = true;
                        break;
                    }
                }
                if (isbroken) {
                    if (excludeBroken) {
                        cerr << var.sequenceName << ":" << var.position << ", sample " << s->first << " has GL > 0" << endl;
                        break;
                    } else {
                        cerr << "VCF record @ " << var.sequenceName << ":" << var.position << ", sample " << s->first << " has GL > 0, not processing, but outputting" << endl;
                        continue;
                    }
                }

                // normalize GLs to -10 min 0 max using division by max and bounding at -10
                double minGL = 0;
                for (vector<double>::iterator g = gls.begin(); g != gls.end(); ++g) {
                    if (*g < minGL) minGL = *g;
                }
                double maxGL = minGL;
                for (vector<double>::iterator g = gls.begin(); g != gls.end(); ++g) {
                    if (*g > maxGL) maxGL = *g;
                }
                // modify gls
                for (vector<double>::iterator g = gls.begin(); g != gls.end(); ++g) {
                    *g = max(glBound, *g - maxGL);
                }

                // and pack back into GL field
                glstrs.clear();
                for (vector<double>::iterator g = gls.begin(); g != gls.end(); ++g) {
                    glstrs.push_back(convert(*g));
                }
            }
        }
        if (excludeBroken && isbroken) {
            cerr << "excluding VCF record @ " << var.sequenceName << ":" << var.position << " due to GLs > 0" << endl;
        } else {
            cout << var << endl;
        }
    }

    return 0;

}

