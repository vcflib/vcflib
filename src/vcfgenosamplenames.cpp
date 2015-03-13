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

    variantFile.addHeaderLine("##FORMAT=<ID=SN,Number=1,Type=String,Description=\"The name of the sample.\">");

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        var.format.push_back("SN");
        for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin();
             s != var.samples.end(); ++s) {
            s->second["SN"].clear();
            s->second["SN"].push_back(s->first);
        }
        cout << var << endl;
    }

    return 0;

}

