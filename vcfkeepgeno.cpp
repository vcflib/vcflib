#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>
#include <set>

using namespace std;
using namespace vcf;

int main(int argc, char** argv) {

    if (argc < 3) {
        cerr << "usage: " << argv[0] << " <vcf file> [FIELD1] [FIELD2] ..." << endl
             << "outputs each record in the vcf file, removing FORMAT fields not listed"
	     << "on the command line from sample specifications in the output"
	     << endl;
        return 1;
    }

    string filename = argv[1];

    set<string> fieldsToKeep;
    for (int i = 2; i < argc; ++i) {
        fieldsToKeep.insert(argv[i]);
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

    // write the header
    cout << variantFile.header;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
	vector<string> fieldsToErase;
	for (map<string, vector<string> >::iterator i = var.info.begin(); i != var.info.end(); ++i) {
	    if (!fieldsToKeep.count(i->first)) {
		fieldsToErase.push_back(i->first);
	    }
	}
	for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
	    map<string, vector<string> >& sample = s->second;
	    for (vector<string>::iterator f = fieldsToErase.begin(); f != fieldsToErase.end(); ++f) {
		sample.erase(*f);
	    }
	}
        cout << var << endl;
    }

    return 0;

}

