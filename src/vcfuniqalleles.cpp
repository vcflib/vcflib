#include "Variant.h"
#include <set>

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

    cout << variantFile.header << endl;

    string lastsn;
    long int lastpos;
    string lastref;
    vector<string> lastalt;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        set<string> alleles;
        vector<string> alleles_to_remove;
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            if (*a != var.ref) {
                if (alleles.find(*a) == alleles.end()) {
                    alleles.insert(*a);
                } else {
                    alleles_to_remove.push_back(*a);
                }
            } else {
                alleles_to_remove.push_back(*a); // same as ref
            }
        }
        for (vector<string>::iterator a = alleles_to_remove.begin(); a != alleles_to_remove.end(); ++a) {
            cerr << "removing " << *a << " from " << var.sequenceName << ":" << var.position << endl;
            var.removeAlt(*a);
        }
        cout << var << endl;
    }

    return 0;

}

