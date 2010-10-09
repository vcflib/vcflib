#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>
#include <getopt.h>

using namespace vcf;


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "options:" << endl 
         << "    -s, --samples    outputs a table of stats / individual in the vcf file (default)" << endl
         << "    -r, --records    outputs a line of statistics over all records in the vcf file" << endl
         << endl;
    exit(0);
}


int main(int argc, char** argv) {

    int c;
    bool statsForSamples = true;
    bool statsForRecords = false;

    if (argc == 1)
        printSummary(argv);

    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"samples",  required_argument, 0, 's'},
            {"records",  required_argument, 0, 'r'},
            //{"length",  no_argument, &printLength, true},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hsr",
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

          case 's':
            statsForSamples = true;
            statsForRecords = false;
            break;
 
          case 'r':
            statsForRecords = true;
            statsForSamples = false;
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

    if (statsForRecords && statsForSamples) {
        cerr << "only one output format may be specified, either -s or -r" << endl;
        exit(1);
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


    map<string, int> sitecount;
    map<string, int> refcount;
    map<string, int> altcount;
    map<string, int> homcount;
    map<string, int> hetcount;
    map<string, int> gqsum;
    map<string, int> dpsum;

    for (vector<string>::iterator s = variantFile.sampleNames.begin(); s != variantFile.sampleNames.end(); ++s) {
        string& sample = *s;
        sitecount[sample] = 0;
        refcount[sample] = 0;
        altcount[sample] = 0;
        homcount[sample] = 0;
        hetcount[sample] = 0;
        gqsum[sample] = 0;
    }

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

        for (map<string, map<string, string> >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {

            string name = s->first;
            map<string, string>& sample = s->second;

            sitecount[name] += 1;

            int gq;
            if (convert(sample["GQ"], gq)) {
                gqsum[name] += gq;
            }

            int dp;
            if (convert(sample["DP"], dp))
                dpsum[name] += dp;

            string& genotype = sample["GT"];
            vector<string> gt = split(genotype, "|/");

            int alt = 0;
            int ref = 0;

            for (vector<string>::iterator g = gt.begin(); g != gt.end(); ++g) {
                if (*g != "0") {
                    ++alt;
                } else {
                    ++ref;
                }
            }

            if (alt != gt.size()) {
                hetcount[name] += alt;
            }

            if (alt == gt.size() || ref == gt.size()) {
                homcount[name] += 1;
            }

            refcount[name] += ref;
            altcount[name] += alt;

        }
    }

    cout << "sample" << "\t"
         << "sitecount" << "\t"
         << "refcount" << "\t"
         << "altcount" << "\t"
         << "homcount" << "\t"
         << "hetcount" << "\t"
         << "avg_gq" << "\t"
         << "avg_dp" << endl;
    for (vector<string>::iterator s = variantFile.sampleNames.begin(); s != variantFile.sampleNames.end(); ++s) {
        string& sample = *s;
        cout << sample << "\t"

             << sitecount[sample] << "\t"
             << refcount[sample] << "\t"
             << altcount[sample] << "\t"
             << homcount[sample] << "\t"
             << hetcount[sample] << "\t"
             << (float) gqsum[sample] / (float) sitecount[sample] << "\t"
             << (float) dpsum[sample] / (float) sitecount[sample]
             << endl;
    }

    return 0;

}

