#include "Variant.h"

using namespace std;
using namespace vcflib;

int main(int argc, char** argv) {

    VariantCallFile variantFile;

    if (argc == 1) {
        cerr << "usage: " << argv[0] << " [key] [vcf_file]" << endl
             << "Sets QUAL from info field tag keyed by [key]." << endl
             << "The VCF file may be omitted and read from stdin." << endl
             << "The average of the field is used if it contains multiple values." << endl;
        return 1;
    }

    string key = argv[1];

    if (argc > 2) {
        string filename = argv[2];
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
        vector<string>& ivs = var.info[key];
        double vs = 0;
        for (vector<string>::iterator i = ivs.begin();
             i != ivs.end(); ++i) {
            double v;
            convert(*i, v);
            vs += v;
        }
        var.quality = vs / (double) ivs.size();
        cout << var << endl;
    }

    return 0;

}

