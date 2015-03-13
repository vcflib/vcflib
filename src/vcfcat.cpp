#include "Variant.h"

using namespace std;
using namespace vcflib;

int main(int argc, char** argv) {

    if (argc == 1) {
        cout << "usage: " << argv[0] << " [file1] [file2] ... [fileN]" << endl
             << "Concatenates VCF files." << endl;
        return 0;
    } else {
        for (int i = 1; i < argc; ++i) {
            VariantCallFile variantFile;
            string filename = argv[i];
            variantFile.open(filename);
            if (!variantFile.is_open()) {
                cerr << "could not open " << argv[i] << endl;
                return 1;
            }
            if (i == 1) {
                cout << variantFile.header << endl;
            }
            Variant var(variantFile);
            while (variantFile.getNextVariant(var)) {
                cout << var << endl;
            }
        }
    }

    return 0;

}

