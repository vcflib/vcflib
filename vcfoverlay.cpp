#include "Variant.h"
#include <getopt.h>

using namespace std;
using namespace vcf;


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [<vcf file> ...]" << endl
         << endl
         << "options:" << endl 
         << "    -h, --help       this dialog" << endl
         << endl
         << "Overlays records in the input vcf files in the order in which they appear." << endl;
    exit(0);
}

int main(int argc, char** argv) {

    if (argc == 1)
        printSummary(argv);

    int c;
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };
        int option_index = 0;

        c = getopt_long (argc, argv, "h",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
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

    // idea here is to shadow-merge
    // records from the VCF files, which are provided in order of desired merge

    map<int, pair<VariantCallFile*, Variant > > variantFiles;
    map<string, map<long int, map<int, string> > > linesByPrecedence;
    int i = optind;

    if (optind < argc - 1) {
	while (i < argc) {
	    int index = i++;
	    VariantCallFile*& variantFile = variantFiles[index].first;
	    Variant& var = variantFiles[index].second;
	    string inputFilename = argv[optind++];
	    variantFile = new VariantCallFile;
	    variantFile->open(inputFilename);
	    var.setVariantCallFile(variantFile);
	    if (!variantFile->is_open()) {
		cout << "could not open VCF file" << endl;
		exit(1);
	    } else {
		if (variantFile->getNextVariant(var)) {
		    linesByPrecedence[var.sequenceName][var.position][index] = variantFile->line;
		}
	    }
	}
    } else {
        cerr << "no input files specified" << endl;
        exit(1);
    }

    cout << variantFiles.begin()->second.first->header << endl;

    while (!linesByPrecedence.empty()) {
        // get the lowest entry in the buffer of observed lines
        // print the first line
        // get the next variant from that file, put it back into the map
        const string& lowestChrom = linesByPrecedence.begin()->first;
        const long int lowestPosition = linesByPrecedence.begin()->second.begin()->first;
        map<int, string>& lowestLine = linesByPrecedence.begin()->second.begin()->second;
        cout << lowestLine.begin()->second << endl;
        
	for (map<int, string>::iterator g = lowestLine.begin(); g != lowestLine.end(); ++g) {
	    int index = g->first;
	    VariantCallFile& variantFile = *variantFiles[index].first;
	    Variant& var = variantFiles[index].second;
	    if (!variantFile.getNextVariant(var)) {
		variantFiles.erase(index);
	    } else {
		linesByPrecedence[var.sequenceName][var.position][index] = variantFile.line;
	    }
	}
	
        linesByPrecedence[lowestChrom].erase(lowestPosition);
        if (linesByPrecedence[lowestChrom].empty()) {
            linesByPrecedence.erase(lowestChrom);
        }
    }

    // flush the rest of the variant records if there are any

    return 0;
}

