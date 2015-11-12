#include "Variant.h"
#include "split.h"
#include "Fasta.h"
#include <getopt.h>
#include <cmath>

using namespace std;
using namespace vcflib;

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <vcf file>" << endl
         << endl
         << "options:" << endl 
         << "    -f, --fasta-reference REF    FASTA reference file to use to obtain primer sequences." << endl
         << "    -n, --number-of-regions N    The number of desired regions." << endl
         << "    -p, --number-of-positions N  The number of positions per region." << endl
         << "    -o, --offset N               Add an offset to region positioning, to avoid boundary" << endl
         << "                                 related artifacts in downstream processing." << endl
         << "    -l, --overlap N              The number of sites to overlap between regions.  Default 0." << endl
         << "    -s, --separator SEQ          Specify string to use to separate region output.  Default '-'" << endl
         << endl
         << "Generates a list of regions, e.g. chr20:10..30 using the variant" << endl
         << "density information provided in the VCF file to ensure that the regions have" << endl
         << "even numbers of variants.  This can be use to reduce the variance in runtime" << endl
         << "when dividing variant detection or genotyping by genomic coordinates." << endl;
    exit(0);
}


struct Region {
    long int start;
    long int end;
    int positions;
    Region() : start(0), end(0), positions(0) { }
    Region(long int s, long int e)
        : start(s), end(e), positions(0) { }
};


int main(int argc, char** argv) {

    int c;
    string fastaRef;
    bool keepFailures = false;
    bool excludeFailures = false;
    int number_of_regions = 1;
    int number_of_positions = 0;
    int offset = 0;
    int overlap = 0;
    string regionSplitSeq = "-";

    if (argc == 1)
        printSummary(argv);

    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"fasta-reference",  required_argument, 0, 'f'},
                {"number-of-regions",  required_argument, 0, 'n'},
                {"number-of-positions",  required_argument, 0, 'p'},
                {"offset",  required_argument, 0, 'o'},
                {"overlap",  required_argument, 0, 'l'},
                {"separator",  required_argument, 0, 's'},
                //{"length",  no_argument, &printLength, true},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hf:n:o:l:s:p:",
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
            fastaRef = optarg;
            break;

        case 'n':
            number_of_regions = atoi(optarg);
            break;

        case 'p':
            number_of_positions = atoi(optarg);
            break;

        case 'o':
            offset = atoi(optarg);
            break;

        case 'l':
            overlap = atoi(optarg);
            break;

        case 's':
            regionSplitSeq = optarg;
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

    if (fastaRef.empty()) {
        cerr << "a FASTA reference sequence must be specified" << endl;
        exit(1);
    }

    FastaReference ref;
    ref.open(fastaRef);

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

    map<string, vector<Region> > positions_by_chrom;
    int total_positions = 0;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        int refstart = var.position - 1; // convert to 0-based
        positions_by_chrom[var.sequenceName].push_back(Region(refstart + offset, refstart + offset + var.ref.size()));
        ++total_positions;
    }

    int positions_per_region;
    if (number_of_positions) {
        if (number_of_positions - overlap < 0) {
            cerr << "overlap is greater than the number of positions per region!" << endl;
            exit(1);
        } else {
            positions_per_region = number_of_positions - overlap;
        }
    } else {
        positions_per_region = ceil((double) total_positions / (double) number_of_regions);
    }
    //cerr << positions_per_region << "=" << total_positions << "/" << number_of_regions << "+" << overlap << endl;

    // todo, update routine to allow overlaps

    for (map<string, vector<Region> >::iterator s = positions_by_chrom.begin();
         s != positions_by_chrom.end(); ++s) {
        //pair<long int, long int> current_region;
        Region current_region;
        for (vector<Region>::iterator p = s->second.begin(); p != s->second.end(); ++p) {
            if (current_region.positions < positions_per_region + overlap) {
                current_region.end = p->end;
                current_region.positions++;
            } else {
                cout << s->first << ":" << current_region.start << regionSplitSeq << current_region.end << endl;
                vector<Region>::iterator l = max(s->second.begin(), p-overlap-1);
                current_region.start = l->end;
                current_region.end = p->end;
                current_region.positions = overlap + 1;
            }
        }
        // get refseq size, use as end coordinate for last region in target
        current_region.end = ref.sequenceLength(s->first);
        cout << s->first << ":" << current_region.start << regionSplitSeq << current_region.end << endl;
    }

    return 0;

}

