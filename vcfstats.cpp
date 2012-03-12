#include "Variant.h"
#include "split.h"
#include "convert.h"
#include <getopt.h>

using namespace std;
using namespace vcf;

bool isTransition(const string& ref, const string& alt) {
    if (((ref == "A" && alt == "G") || (ref == "G" && alt == "A")) ||
        ((ref == "C" && alt == "T") || (ref == "T" && alt == "C"))) {
        return true;
    } else {
        return false;
    }
}

bool isDeamination(const string& ref, const string& alt) {
    if ((ref == "G" && alt == "A") ||
        (ref == "C" && alt == "T")) {
        return true;
    } else {
        return false;
    }
}

bool isAmination(const string& ref, const string& alt) {
    if ((ref == "A" && alt == "G") ||
        (ref == "T" && alt == "C")) {
        return true;
    } else {
        return false;
    }
}

class AlleleStats {
public:
    int transitions;
    int transversions;
    int deaminations;
    int aminations;
    int mismatches;
    int insertions;
    int deletions;
    //AlleleStats(int ts, int tv, int da, int am, int mm)
    AlleleStats(void)
        : transitions(0)
        , transversions(0)
        , deaminations(0)
        , aminations(0)
        , mismatches(0)
        , insertions(0)
        , deletions(0)
    { }
};

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "    -r, --region          specify a region on which to target the stats, requires a BGZF" << endl
         << "                          compressed file which has been indexed with tabix.  any number of" << endl
         << "                          regions may be specified." << endl
         << "    -a, --add-info        add the statistics intermediate information to the VCF file," << endl
         << "                          writing out VCF records instead of summary statistics" << endl
	 << "    -l, --no-length-frequency   don't out the indel and mnp length-frequency spectra" << endl
         << endl
         << "Prints statistics about variants in the input VCF file." << endl;
}


int main(int argc, char** argv) {

    vector<string> regions;
    bool addTags = false;
    bool lengthFrequency = true;

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"region", required_argument, 0, 'r'},
            {"add", no_argument, 0, 'a'},
            {"no-length-frequency", no_argument, 0, 'l'},
            //{"length",  no_argument, &printLength, true},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hlar:",
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
		
	    case 'l':
		lengthFrequency = false;
		break;
		
	    case 'a':
		addTags = true;
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

    if (addTags) {
        variantFile.addHeaderLine("##INFO=<ID=transitions,Number=A,Type=Integer,Description=\"Total number of transitions in the alternate allele\">");
        variantFile.addHeaderLine("##INFO=<ID=transversions,Number=A,Type=Integer,Description=\"Total number of transversions in the alternate allele\">");
        variantFile.addHeaderLine("##INFO=<ID=deaminations,Number=A,Type=Integer,Description=\"Total number of deaminations in the alternate allele\">");
        variantFile.addHeaderLine("##INFO=<ID=aminations,Number=A,Type=Integer,Description=\"Total number of aminations in the alternate allele\">");
        variantFile.addHeaderLine("##INFO=<ID=mismatches,Number=A,Type=Integer,Description=\"Total number of mismatches in the alternate allele\">");
        variantFile.addHeaderLine("##INFO=<ID=insertions,Number=A,Type=Integer,Description=\"Total number of inserted bases in the alternate allele\">");
        variantFile.addHeaderLine("##INFO=<ID=deletions,Number=A,Type=Integer,Description=\"Total number of deleted bases in the alternate allele\">");
        cout << variantFile.header << endl;
    }

    Variant var(variantFile);

    vector<string>::iterator regionItr = regions.begin();

    int variantAlleles = 0;
    int uniqueVariantAlleles = 0;
    int variantSites = 0;
    int snps = 0;
    int transitions = 0;
    int transversions = 0;
    int deaminations = 0;
    int aminations = 0;
    int totalinsertions = 0;
    int totaldeletions = 0;
    int insertedbases = 0;
    int deletedbases = 0;
    int totalmnps = 0;
    int totalcomplex = 0;
    int mismatchbases = 0;
    int mnpbases = 0;
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
            map<VariantAllele, vector<string> > uniqueVariants;
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                string& alternate = *a;
                vector<VariantAllele>& vav = alternates[alternate];
                if (vav.size() > 1) {
                    ++totalcomplex;
                }
                for (vector<VariantAllele>::iterator v = vav.begin(); v != vav.end(); ++v) {
                    uniqueVariants[*v].push_back(alternate);
                }
            }
            uniqueVariantAlleles += uniqueVariants.size();
            variantAlleles += var.alt.size();
            map<string, AlleleStats> alleleStats;

            for (map<VariantAllele, vector<string> >::iterator v = uniqueVariants.begin(); v != uniqueVariants.end(); ++v) {
                const VariantAllele& va = v->first;
                vector<string>& alternates = v->second;
                if (!addTags) { // don't add any tag information if we're not going to output it
                    alternates.clear();
                }

                if (va.ref.size() == va.alt.size()) {
                    if (va.ref.size() == 1) {
                        ++snps;
                        ++mismatchbases;
                        for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                            ++alleleStats[*a].mismatches;
                        }
                        if (isTransition(va.ref, va.alt)) {
                            ++transitions;
                            for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                                ++alleleStats[*a].transitions;
                            }
                        } else {
                            ++transversions;
                            for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                                ++alleleStats[*a].transversions;
                            }
                        }
                        if (isAmination(va.ref, va.alt)) {
                            ++aminations;
                            for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                                ++alleleStats[*a].aminations;
                            }
                        }
                        if (isDeamination(va.ref, va.alt)) {
                            ++deaminations;
                            for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                                ++alleleStats[*a].deaminations;
                            }
                        }
                    } else {
                        ++totalmnps;
                        ++mnps[va.alt.size()]; // not entirely correct
                        for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                            alleleStats[*a].mismatches += va.alt.size();
                        }
                        string::const_iterator r = va.ref.begin();
                        for (string::const_iterator a = va.alt.begin(); a != va.alt.end(); ++a, ++r) {
                            string rstr = string(1, *r);
                            string astr = string(1, *a);
                            if (rstr == astr) {
                                continue;
                            }
                            if (isTransition(rstr, astr)) {
                                ++transitions;
                                for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                                    ++alleleStats[*a].transitions;
                                }
                            } else {
                                ++transversions;
                                for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                                    ++alleleStats[*a].transversions;
                                }
                            }
                            if (isAmination(rstr, astr)) {
                                ++aminations;
                                for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                                    ++alleleStats[*a].aminations;
                                }
                            }
                            if (isDeamination(rstr, astr)) {
                                ++deaminations;
                                for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                                    ++alleleStats[*a].deaminations;
                                }
                            }
                            ++mismatchbases;
                            ++mnpbases;
                        }
                    }
                } else if (va.ref.size() > va.alt.size()) {
                    int diff = va.ref.size() - va.alt.size();
                    deletedbases += diff;
                    ++totaldeletions;
                    ++deletions[diff];
                    for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                        alleleStats[*a].deletions += diff;
                    }
                } else {
                    int diff = va.alt.size() - va.ref.size();
                    insertedbases += diff;
                    ++totalinsertions;
                    ++insertions[diff];
                    for (vector<string>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                        alleleStats[*a].insertions += diff;
                    }
                }
            }
            if (addTags) {
                for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                    var.info["mismatches"].push_back(convert(alleleStats[*a].mismatches));
                    var.info["insertions"].push_back(convert(alleleStats[*a].insertions));
                    var.info["deletions"].push_back(convert(alleleStats[*a].deletions));
                    var.info["transitions"].push_back(convert(alleleStats[*a].transitions));
                    var.info["transversions"].push_back(convert(alleleStats[*a].transversions));
                    var.info["deaminations"].push_back(convert(alleleStats[*a].deaminations));
                    var.info["aminations"].push_back(convert(alleleStats[*a].aminations));
                }
                cout << var << endl;
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

    if (!addTags) {
        cout << "total variant sites:\t" << variantSites << endl
             << "total variant alleles:\t" << variantAlleles << endl
             << "unique variant alleles:\t" << uniqueVariantAlleles << endl
             << endl
             << "snps:\t" << snps << endl
             << "mnps:\t" << totalmnps << endl
             << "indels:\t" << totalinsertions + totaldeletions << endl
             << "complex:\t" << totalcomplex << endl
             << endl
             << "mismatches:\t" << mismatchbases << endl
             << endl
             << "ts/tv ratio:\t" << (double) transitions / (double) transversions << endl
             << "deamination ratio:\t" << (double) deaminations / aminations << endl;
	if (lengthFrequency) {
	    cout << endl
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
	}
        cout << endl
             << "insertion alleles / deletion alleles:\t" << (double) totalinsertions / (double) totaldeletions << endl
             << "inserted bases / deleted bases:\t" << (double) insertedbases / (double) deletedbases << endl
             << endl;
	if (lengthFrequency) {
	    cout << "mnp length frequency distribution" << endl
		 << "length\tcount" << endl;
	    for (int i = 2; i <= maxmnp; ++i) {
		int mnp = mnps[i];
		cout << i << "\t"
		     << (mnp > 0 ? convert(mnp) : "")
		     << endl;
	    }
	}
        cout << "total bases in mnps:\t" << mnpbases << endl;

        /*
        cout << "complex event frequency distribution" << endl
             << "length\tcount" << endl;
        for (map<int, int>::iterator i = complexsubs.begin(); i != complexsubs.end(); ++i) {
            cout << i->first << "\t" << i->second << endl;
        }
        */
    }

    return 0;

}

