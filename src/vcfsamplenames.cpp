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

    for (vector<string>::iterator sample = variantFile.sampleNames.begin();
            sample != variantFile.sampleNames.end(); ++sample) {
        cout << *sample << endl;
    }

    return 0;

}

