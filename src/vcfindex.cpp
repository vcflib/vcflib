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

    string idname = "id";
    long int uid = 0;

    variantFile.addHeaderLine("##INFO=<ID="+idname+",Number=A,Type=Integer,Description=\"Unique numerical identifier of allele in file.\">");
    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        vector<string>& idxs = var.info[idname];
        idxs.clear();
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            idxs.push_back(convert(uid++));
        }
        cout << var << endl;
    }

    return 0;

}

