#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

    if (argc > 1) {
        cerr << "usage: " << argv[0] << " <[vcf file]" << endl
             << "modifies the genotypes field to provide the literal alleles rather than indexes" << endl;
        return 1;
    }

    VariantCallFile variantFile;

    variantFile.open(std::cin);

    if (!variantFile.is_open()) {
        return 1;
    }

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        map<string, map<string, vector<string> > >::iterator s     = var.samples.begin(); 
        map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();
        
        for (; s != sEnd; ++s) {
            map<string, vector<string> >& sample = s->second;
            vector<string>& gtstrs = sample["GT"];
            string& genotype = gtstrs.front();
            vector<string> gt = split(genotype, "|/");
            
            // report the sample and it's genotype
            stringstream o;
            for (vector<string>::iterator g = gt.begin(); g != gt.end(); ++g) {
                int index = atoi(g->c_str());
                o << var.alleles[index];
                if (g != (gt.end()-1)) o << "/";
            }
            gtstrs.clear();
            gtstrs.push_back(o.str());
        }
        cout << var << endl;
    }
    return 0;

}

