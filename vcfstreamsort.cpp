#include "Variant.h"

using namespace std;
using namespace vcf;

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

    cout << variantFile.header << endl;

    map<long int, vector<Variant> > records;
    long int back = 0;
    int sortSitesWindow = 100;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
	records[var.position].push_back(var);
	if (records.size() > sortSitesWindow) {
	    vector<Variant>& vars = records.begin()->second;
	    for (vector<Variant>::iterator v = vars.begin(); v != vars.end(); ++v) {
		cout << *v << endl;
	    }
	    records.erase(records.begin());
	}
    }
    for (map<long int, vector<Variant> >::iterator r = records.begin(); r != records.end(); ++r) {
	for (vector<Variant>::iterator v = r->second.begin(); v != r->second.end(); ++v) {
	    cout << *v << endl;
	}
    }

    return 0;

}

