#include "Variant.h"
#include "BedReader.h"
#include "intervaltree/IntervalTree.h"
#include <getopt.h>
#include "fastahack/Fasta.h"
#include <algorithm>
#include <list>
#include <set>

using namespace std;
using namespace vcf;


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [<vcf file>]" << endl
         << endl
         << "options:" << endl 
	 << "    -w, --ref-window-size N      align using this many bases flanking each side of the reference allele" << endl
	 << "    -s, --alt-window-size N      align using this many flanking bases from the reference around each alternate allele" << endl
	 << "    -r, --reference FILE         FASTA reference file, required with -i and -u" << endl
	 << "    -m, --match-score N          match score for SW algorithm" << endl
	 << "    -x, --mismatch-score N       mismatch score for SW algorithm" << endl
	 << "    -o, --gap-open-penalty N     gap open penalty for SW algorithm" << endl
	 << "    -e, --gap-extend-penalty N   gap extension penalty for SW algorithm" << endl
	 << "    -z, --entropy-gap-open       use entropy scaling for the gap open penalty" << endl
	 << "    -R, --repeat-gap-extend N    penalize non-repeat-unit gaps in repeat sequence" << endl
         << endl
	 << "For each alternate allele, attempt to realign against the reference with lowered gap open penalty." << endl
	 << "If realignment is possible, adjust the cigar and reference/alternate alleles." << endl;
    exit(0);
}

int main(int argc, char** argv) {

    string vcfFileName;
    string fastaFileName;
    int windowsize = 100;
    bool includePreviousBaseForIndels = false;
    bool useMNPs = true;
    int altwindowsize = 50;

    // constants for SmithWaterman algorithm
    float matchScore = 10.0f;
    float mismatchScore = -9.0f;
    float gapOpenPenalty = 15.0f;
    float gapExtendPenalty = 6.66f;

    bool useEntropy = false;
    bool useRepeatGapExtendPenalty = false;
    float repeatGapExtendPenalty = 1;

    if (argc == 1)
        printSummary(argv);

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
	    {"ref-window-size", required_argument, 0, 'w'},
	    {"reference", required_argument, 0, 'r'},
	    {"match-score", required_argument, 0, 'm'},
	    {"mismatch-score", required_argument, 0, 'x'},
	    {"gap-open-penalty", required_argument, 0, 'o'},
	    {"gap-extend-penalty", required_argument, 0, 'e'},
	    {"alt-window-size", required_argument, 0, 's'},
	    {"entropy-gap-open", no_argument, 0, 'z'},
	    {"repeat-gap-extend", no_argument, 0, 'R'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hzw:r:m:x:o:e:s:R:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

	    case 'w':
		windowsize = atoi(optarg);
		break;

	    case 'r':
		fastaFileName = string(optarg);
		break;

            case 'h':
                printSummary(argv);
                break;

	    case 'm':
		matchScore = atof(optarg);
	        break;

	    case 'x':
		mismatchScore = atof(optarg);
	        break;

	    case 'o':
		gapOpenPenalty = atof(optarg);
	        break;

	    case 'e':
		gapExtendPenalty = atof(optarg);
	        break;

	    case 's':
		altwindowsize = atoi(optarg);
		break;

	    case 'z':
		useEntropy = true;
		break;

	    case 'R':
		useRepeatGapExtendPenalty = true;
		repeatGapExtendPenalty = atof(optarg);
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
        cerr << "could not open VCF file" << endl;
        exit(1);
    }

    FastaReference freference;
    if (fastaFileName.empty()) {
	cerr << "a reference is required" << endl;
	exit(1);
    } else {
	freference.open(fastaFileName);
    }

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
	cout << endl;
	cout << var << endl;
	map<string, vector<VariantAllele> > variantAlleles;
	for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
	    cout << endl;
	    //for (int i = 10; i < 100; ++i) {

	    //altwindowsize = i;

	    // try to remap locally

	    string reference = freference.getSubSequence(var.sequenceName, var.position - 1 - windowsize, windowsize * 2 + var.ref.size());
	    
	    // passed to sw align
	    unsigned int referencePos;
	    string cigar;

	    string& alternate = *a;

	    vector<VariantAllele>& variants = variantAlleles[alternate];

	    string alternateQuery = reference.substr(windowsize - altwindowsize, altwindowsize) + alternate + reference.substr(reference.size() - windowsize, altwindowsize);

	    //cout << "REF:\t" << reference << endl;
	    //cout << "ALT:\t" << string(windowsize - altwindowsize, ' ') << alternateQuery << endl;
	    
	    CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
	    if (useEntropy) sw.EnableEntropyGapPenalty(1);
	    if (useRepeatGapExtendPenalty) sw.EnableRepeatGapExtensionPenalty(repeatGapExtendPenalty);
	    sw.Align(referencePos, cigar, reference, alternateQuery);

	    // left-realign the alignment...

	    int altpos = 0;
	    int refpos = 0;
	    int len;
	    string slen;
	    vector<pair<int, char> > cigarData;

	    string ref = reference.substr(referencePos);

	    stringstream refss;
	    stringstream altss;

	    cout << cigar << endl;
	    for (string::iterator c = cigar.begin(); c != cigar.end(); ++c) {
		switch (*c) {
                case 'I':
                    len = atoi(slen.c_str());
                    slen.clear();
                    cigarData.push_back(make_pair(len, *c));
		    altss << alternateQuery.substr(altpos, len);
		    refss << string(len, '-');
                    altpos += len;
		    break;
		case 'D':
                    len = atoi(slen.c_str());
                    slen.clear();
                    cigarData.push_back(make_pair(len, *c));
		    refss << ref.substr(refpos, len);
		    altss << string(len, '-');
                    refpos += len;
		    break;
		case 'M':
		    len = atoi(slen.c_str());
		    slen.clear();
		    cigarData.push_back(make_pair(len, *c));
		    refss << ref.substr(refpos, len);
		    altss << alternateQuery.substr(altpos, len);
		    refpos += len;
		    altpos += len;
		    break;
		case 'S':
                    len = atoi(slen.c_str());
                    slen.clear();
                    cigarData.push_back(make_pair(len, *c));
		    refss << ref.substr(refpos, len);
 		    //altss << alternateQuery.substr(altpos, len); // TODO deal with soft clipping, weird behavior
                    refpos += len;
                    altpos += len;
		    break;
		default:
                    len = 0;
                    slen += *c;
		    break;
		}
	    }

	    cout << "ref:\t" << refss.str() << endl;
	    cout << "alt:\t" << altss.str() << endl;

//	    }
	}

	
    }

    return 0;

}

