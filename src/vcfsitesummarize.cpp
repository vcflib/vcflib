#include "Variant.h"

using namespace std;
using namespace vcflib;

int main(int argc, char** argv) {

    VariantCallFile variantFile;

    if (argc > 1) {
        string filename = argv[1];
        variantFile.open(filename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    // obtain all possible field names
    vector<string> infofields;
    vector<string> infoflags;

    for (map<string, VariantFieldType>::iterator i = variantFile.infoTypes.begin(); i != variantFile.infoTypes.end(); ++i) {
        if (variantFile.infoCounts[i->first] != ALLELE_NUMBER) {
            if (i->second == FIELD_BOOL) {
                infoflags.push_back(i->first);
            } else {
                infofields.push_back(i->first);
            }
        }
    }

    // write header

    // defaults
    cout << "CHROM\tPOS\tID\tREF\tQUAL\tFILTER";
    
    // configurable info field
    for (vector<string>::iterator i = infofields.begin(); i != infofields.end(); ++i) {
        cout << "\t" << *i;
    }
    for (vector<string>::iterator i = infoflags.begin(); i != infoflags.end(); ++i) {
        cout << "\t" << *i;
    }
    cout << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

	cout << var.sequenceName << "\t"
	     << var.position << "\t"
	     << var.id << "\t"
	     << var.ref << "\t"
	     << var.quality << "\t"
	     << var.filter;

	for (vector<string>::iterator i = infofields.begin(); i != infofields.end(); ++i) {
	    vector<string> value;
	    string& name = *i;
	    map<string, vector<string> >::iterator f = var.info.find(name);
	    if (f != var.info.end()) {
            value = f->second;
            if (value.size() == 1) {
                cout << "\t" << value.front();
            } else {
                cout << "\t"; // null
            }
	    } else {
            cout << "\t"; // null
	    }
	}

	for (vector<string>::iterator i = infoflags.begin(); i != infoflags.end(); ++i) {
	    string value;
	    string& name = *i;
	    map<string, bool>::iterator f = var.infoFlags.find(name);
	    cout << "\t";
	    if (f != var.infoFlags.end()) {
            cout << 1;
	    } else {
            cout << 0;
	    }
	}
	
	cout << endl;
	
    }
    
    return 0;

}

