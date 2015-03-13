#include "Variant.h"
#include "convert.h"
#include <vector>

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

    variantFile.addHeaderLine("##INFO=<ID=length,Number=A,Type=Integer,Description=\"length(ALT) - length(REF) for each ALT\">");
    variantFile.addHeaderLine("##INFO=<ID=length.ref,Number=1,Type=Integer,Description=\"length(REF)\">");
    variantFile.addHeaderLine("##INFO=<ID=length.alt,Number=A,Type=Integer,Description=\"length(ALT) for each ALT\">");
    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        vector<string>& lengths = var.info["length"];
        lengths.clear();
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            lengths.push_back(convert((int) a->size() - (int) var.ref.size()));
        }
        vector<string>& lengthsRef = var.info["length.ref"];
        lengthsRef.clear();
        lengthsRef.push_back(convert(var.ref.size()));
        vector<string>& lengthsAlt = var.info["length.alt"];
        lengthsAlt.clear();
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            lengthsAlt.push_back(convert((int) a->size()));
        }
        cout << var << endl;
    }

    return 0;

}

