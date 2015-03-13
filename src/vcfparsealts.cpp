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

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        map<string, vector<VariantAllele> > variants = var.parsedAlternates();
	cout << var << endl;
        for (map<string, vector<VariantAllele> >::iterator va = variants.begin(); va != variants.end(); ++va) {
            cout << " ( " << va->first << " :: ";
            vector<VariantAllele>& vars = va->second;
            vector<VariantAllele>::iterator g = vars.begin();
            for (; g != vars.end(); ++g) {
                cout << *g << "; ";
            }
            cout << " ) ";
        }
        cout << endl;
    }

    return 0;

}

