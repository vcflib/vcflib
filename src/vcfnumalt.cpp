#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
using namespace vcflib;

int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "outputs a VCF stream where NUMALT has been generated for each record using sample genotypes" << endl;
        return 1;
    }

    string filename = argv[1];

    VariantCallFile variantFile;
    if (filename == "-") {
        variantFile.open(std::cin);
    } else {
        variantFile.open(filename);
    }

    if (!variantFile.is_open()) {
        cerr << "could not open " << filename << endl;
        return 1;
    }

    Variant var(variantFile);

    // remove header lines we're going to add
    variantFile.removeInfoHeaderLine("NUMALT");

    // and add them back, so as not to duplicate them if they are already there
    variantFile.addHeaderLine("##INFO=<ID=NUMALT,Number=1,Type=Integer,Description=\"Total number of segregating alternate alleles at the loci\">");

    // write the new header
    cout << variantFile.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        stringstream na;
        na << var.alt.size();
        var.info["NUMALT"].clear();
        var.info["NUMALT"].push_back(na.str());
        cout << var << endl;
    }

    return 0;

}

