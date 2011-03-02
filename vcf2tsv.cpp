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
    // obtain all possible field names
    vector<string> infofields;
    vector<string> infoflags;

    for (map<string, string>::iterator i = variantFile.infoTypes.begin(); i != variantFile.infoTypes.end(); ++i) {
        if (i->second == "Flag") {
            infoflags.push_back(i->first);
        } else {
            infofields.push_back(i->first);
        }
    }

    // write header

    // defaults
    cout << "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t";
    
    // configurable info field
    for (vector<string>::iterator i = infofields.begin(); i != infofields.end(); ++i) {
        if (i != infofields.begin()) {
            cout << "\t";
        }
        cout << *i;
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
             << var.ref << "\t";
        var.printAlt(cout);
        cout << "\t"
             << var.quality << "\t"
             << var.filter;

        for (vector<string>::iterator i = infofields.begin(); i != infofields.end(); ++i) {
            string value;
            string& name = *i;
            map<string, string>::iterator f = var.info.find(name);
            if (f != var.info.end()) {
                value = f->second;
            }
            cout << "\t" << value;
        }

        for (vector<string>::iterator i = infoflags.begin(); i != infoflags.end(); ++i) {
            string value;
            string& name = *i;
            map<string, bool>::iterator f = var.infoFlags.find(name);
            cout << "\t";
            if (f != var.infoFlags.end()) {
                cout << "TRUE";
            } else {
                cout << "FALSE";
            }
        }

        cout << endl;

    }

    return 0;

}

