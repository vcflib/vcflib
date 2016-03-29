#include "Variant.h"
#include <getopt.h>
#include "gpatInfo.hpp"


using namespace std;
using namespace vcflib;


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [<vcf file> ...]" << endl
         << endl
         << "options:" << endl 
         << "    -h, --help       this dialog" << endl
	 << "    -v, --version    prints version" << endl
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
	  {"version", no_argument, 0, 'v'},
	  {0, 0, 0, 0}
	};
      int option_index = 0;
      
      c = getopt_long (argc, argv, "hv",
		       long_options, &option_index);
      
      if (c == -1){
	break;
      }
      switch (c) {
      case 'h':
	{
	  printSummary(argv);
	  break;
	}
      case 'v':
	{
	  printBasicVersion();
	  exit(0);
	} 
      case '?':
	{
	  printSummary(argv);
	  exit(1);
	  break;
	}
      default:
	abort ();
      }
    }

    // idea here is to shadow-merge
    // records from the VCF files, which are provided in order of desired merge

    map<int, pair<VariantCallFile*, Variant > > variantFiles;
    map<string, map<long int, map<string, map<int, string> > > > linesByPrecedence;
    int i = optind;

    if (!(optind < argc - 1)) {
        cerr << "more than one input file must be specified" << endl;
        exit(1);
    }

	while (i < argc) {
	    int index = i++;
	    VariantCallFile*& variantFile = variantFiles[index].first;
	    Variant& var = variantFiles[index].second;
	    string inputFilename = argv[optind++];
	    variantFile = new VariantCallFile;
        try {
            if (!variantFile->open(inputFilename)) {
                cerr << "vcfoverlay could not open VCF file " << inputFilename << endl;
                --index;
            } else {
                var.setVariantCallFile(variantFile);
                while (variantFile->getNextVariant(var)) {
                    linesByPrecedence[var.sequenceName][var.position][var.vrepr()][index] = variantFile->line;
                }
            }
        } catch (...) {
            cerr << "vcfoverlay encountered errors when opening " << inputFilename << endl;
        }
    }
    
    cout << variantFiles.begin()->second.first->header << endl;

    while (!linesByPrecedence.empty()) {
        // get the lowest entry in the buffer of observed lines
        // print the first line
        // get the next variant from that file, put it back into the map
        const string& lowestChrom = linesByPrecedence.begin()->first;
        const long int lowestPosition = linesByPrecedence.begin()->second.begin()->first;
        map<string, map<int, string> >& pos = linesByPrecedence.begin()->second.begin()->second;
        for (map<string, map<int, string> >::iterator m = pos.begin(); m != pos.end(); ++m) {
            cout << m->second.begin()->second << endl;
        }
        linesByPrecedence[lowestChrom].erase(lowestPosition);
        
        if (linesByPrecedence[lowestChrom].empty()) {
            linesByPrecedence.erase(lowestChrom);
        }
    }

    // flush the rest of the variant records if there are any

    return 0;
}

