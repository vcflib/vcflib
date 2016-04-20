#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>
#include <set>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

    if (argc < 3) {
        cerr << "usage: " << argv[0] << " <vcf file> [FIELD1] [FIELD2] ..." << endl
             << "outputs each record in the vcf file, removing FORMAT fields not listed "
	     << "on the command line from sample specifications in the output"
	     << endl;
        return 1;
    }

    string filename = argv[1];

    vector<string> newFormat;
    set<string> fieldsToKeep;
    for (int i = 2; i < argc; ++i) {
        fieldsToKeep.insert(argv[i]);
        newFormat.push_back(argv[i]);
    }

    VariantCallFile variantFile;
    if (filename == "-") {
        variantFile.open(std::cin);
    } else {
        variantFile.open(filename);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    Variant var(variantFile);

    vector<string> formatIds = variantFile.formatIds();
    for (vector<string>::iterator i = formatIds.begin(); i != formatIds.end(); ++i) {
        if (!fieldsToKeep.count(*i)) {
            variantFile.removeGenoHeaderLine(*i);
        }
    }

    // write the header
    cout << variantFile.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        var.format = newFormat;
        cout << var << endl;
    }

    return 0;

}

