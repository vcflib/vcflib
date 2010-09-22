#include "Variant.h"


int main(int argc, char** argv) {

    if (argc != 2) {
        return 1;
    }

    string filename = argv[1];

    VariantCallFile variantFile;
    if (!variantFile.openVCF(filename)) {
        return 1;
    }

    cout << variantFile.header << endl;

    Variant var(variantFile.sampleNames);
    while (variantFile.getNextVariant(var)) {
        cout << var << endl;
    }

    return 0;

}

