#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcf;


int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "report the genotypes for each sample, for each variant in the vcf file" << endl;
        return 1;
    }

    string filename = argv[1];

    VariantCallFile variantFile(filename);
    if (!variantFile.is_open()) {
        return 1;
    }

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        map<string, map<string, string> >::iterator s     = var.samples.begin(); 
        map<string, map<string, string> >::iterator sEnd  = var.samples.end();
        for (; s != sEnd; ++s) {
            map<string, string>& sample = s->second;
            string& genotype = sample["GT"];
            vector<string> gt = split(genotype, "|/");
            cout << s->first << ":";
            for (vector<string>::iterator g = gt.begin(); g != gt.end(); ++g) {
                int index = atoi(g->c_str());
                cout << var.alleles[index] << ",";
            }
            cout << "\t";
        }
        cout << "\t" << endl;
    }
    return 0;

}

