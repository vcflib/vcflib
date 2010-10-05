#include "Variant.h"
#include "Split.h"
#include <string>
#include <iostream>

using namespace std;


int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "count the number of alternate alleles in all records in the vcf file" << endl;
        return 1;
    }

    string filename = argv[1];

    VariantCallFile variantFile;
    if (!variantFile.openVCF(filename)) {
        return 1;
    }

    unsigned int alternateAlleleCount = 0;

    Variant var(variantFile.sampleNames, variantFile);
    while (variantFile.getNextVariant(var)) {
        //cout << var << endl;
        for (map<string, map<string, string> >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
            //string& name = s->first;
            map<string, string>& sample = s->second;
            string& genotype = sample["GT"];
            vector<string> gt = split(genotype, "|/");
            int alt = 0;
            for (vector<string>::iterator g = gt.begin(); g != gt.end(); ++g) {
                if (*g != "0")
                    ++alt;
            }
            alternateAlleleCount += alt;
        }
    }

    cout << alternateAlleleCount << endl;

    return 0;

}

