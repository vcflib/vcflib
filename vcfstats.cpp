#include "Variant.h"
#include "split.h"
#include "convert.h"
#include <getopt.h>

using namespace std;
using namespace vcf;

bool isTransition(string& ref, string& alt) {
    if (((ref == "A" && alt == "G") || (ref == "G" && alt == "A")) ||
        ((ref == "C" && alt == "T") || (ref == "T" && alt == "C"))) {
        return true;
    } else {
        return false;
    }
}

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "    -r, --region          specify a region on which to target the stats, requires a BGZF" << endl
         << "                          compressed file which has been indexed with tabix.  any number of" << endl
         << "                          regions may be specified." << endl
         << endl
         << "Prints statistics about variants in the input VCF file." << endl;
}


int main(int argc, char** argv) {

    vector<string> regions;

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"region", required_argument, 0, 'r'},
            //{"length",  no_argument, &printLength, true},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hr:",
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

          case 'h':
            printSummary(argv);
            exit(0);
            break;

          case 'r':
            regions.push_back(optarg);
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
        return 1;
    }

    Variant var(variantFile);

    vector<string>::iterator regionItr = regions.begin();

    int variantAlleles = 0;
    int variantSites = 0;
    int snps = 0;
    int transitions = 0;
    int transversions = 0;
    int totalinsertions = 0;
    int totaldeletions = 0;
    int insertedbases = 0;
    int deletedbases = 0;
    int totalmnps = 0;
    int totalcomplex = 0;
    map<int, int> insertions;
    map<int, int> deletions;
    map<int, int> mnps;
    map<int, int> complexsubs;

    do {

        if (!inputFilename.empty() && !regions.empty()) {
            string regionStr = *regionItr++;
            variantFile.setRegion(regionStr);
        }

        while (variantFile.getNextVariant(var)) {
            ++variantSites;
            map<string, vector<VariantAllele> > alternates = var.parsedAlternates();
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                ++variantAlleles;
                string& alternate = *a;
                if (var.ref.size() == alternate.size()) {
                    if (var.ref.size() == 1) {
                        ++snps;
                        if (isTransition(var.ref, alternate)) {
                            ++transitions;
                        } else {
                            ++transversions;
                        }
                    } else {
                        ++totalmnps;
                        if (alternates[alternate].size() > 1) {
                        } else {
                            VariantAllele& va = alternates[alternate].front();
                            ++mnps[va.alt.size()]; // not entirely correct
                        }
                    }
                } else if (var.ref.size() > alternate.size()) {
                    int diff = var.ref.size() - alternate.size();
                    deletedbases += diff;
                    if (alternates[alternate].size() > 1) {
                        ++totalcomplex;
                        ++complexsubs[-diff];
                    } else {
                        ++totaldeletions;
                        ++deletions[diff];
                    }
                } else {
                    int diff = alternate.size() - var.ref.size();
                    insertedbases += diff;
                    if (alternates[alternate].size() > 1) {
                        ++totalcomplex;
                        ++complexsubs[diff];
                    } else {
                        ++totalinsertions;
                        ++insertions[diff];
                    }
                }
            }
        }

    } while (regionItr != regions.end());

    // find the maximum indel size
    int maxindel = 0;
    for (map<int, int>::iterator i = insertions.begin(); i != insertions.end(); ++i) {
        if (i->first > maxindel) {
            maxindel = i->first;
        }
    }
    for (map<int, int>::iterator i = deletions.begin(); i != deletions.end(); ++i) {
        if (i->first > maxindel) {
            maxindel = i->first;
        }
    }

    // and maximum mnp
    int maxmnp = 0;
    for (map<int, int>::iterator i = mnps.begin(); i != mnps.end(); ++i) {
        if (i->first > maxmnp) {
            maxmnp = i->first;
        }
    }

    // now print the results

    cout << "total variant sites:\t" << variantSites << endl
         << "total variant alleles:\t" << variantAlleles << endl
         << endl
         << "snps:\t" << snps << endl
         << "indels:\t" << totalinsertions + totaldeletions << endl
         << "mnps:\t" << totalmnps << endl
         << "complex:\t" << totalcomplex << endl
         << endl
         << "ts/tv ratio:\t" << (double) transitions / (double) transversions << endl
         << endl
         << "ins/del length frequency distribution" << endl
         << "length\tins\tdel\tins/del" << endl;
    for (int i = 1; i <= maxindel; ++i) {
        int ins = insertions[i];
        int del = deletions[i];
        cout << i << "\t"
             << (ins > 0 ? convert(ins) : "" ) << "\t"
             << (del > 0 ? convert(del) : "") << "\t"
             << (ins > 0 && del > 0 ? convert((double) ins / (double) del) : "")
             << endl;
    }
    cout << endl
         << "insertion alleles / deletion alleles:\t" << (double) totalinsertions / (double) totaldeletions << endl
         << "inserted bases / deleted bases:\t" << (double) insertedbases / (double) deletedbases << endl
         << endl
         << "mnp length frequency distribution" << endl
         << "length\tcount" << endl;
    for (int i = 2; i <= maxmnp; ++i) {
        int mnp = mnps[i];
        cout << i << "\t"
             << (mnp > 0 ? convert(mnp) : "")
             << endl;
    }
    cout << endl;

    cout << "complex event frequency distribution" << endl
         << "length\tcount" << endl;
    for (map<int, int>::iterator i = complexsubs.begin(); i != complexsubs.end(); ++i) {
        cout << i->first << "\t" << i->second << endl;
    }

    return 0;

}

