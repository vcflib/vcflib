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

    int uniqueAlleles = 0;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        uniqueAlleles += var.alleles.size();
    }

    cout << uniqueAlleles << endl;

    return 0;

}

