#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "replaces the null parts of genotypes with the reference alleles" << endl;
        return 1;
    }

    string filename = argv[1];

    VariantCallFile variantFile;

    if (filename == "-") {
        variantFile.open(std::cin);
    } else {
        variantFile.open(filename);
    }

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
            string& genotype = sample["GT"].front();
            string splitter = "/";
            if (genotype.find("/") == string::npos
                || genotype.find("|") == string::npos) {
                // assume haploid
            } else {
                splitter = "|";
            }
            vector<string> gt = split(genotype, "|/");
            // now let's flatten the genotype
            stringstream gts;
            for (vector<string>::iterator g = gt.begin(); g != gt.end(); ++g) {
                if (g != gt.begin()) {
                    gts << splitter;
                }
                if (*g == ".") {
                    gts << 0;
                } else {
                    gts << *g;
                }
            }
            genotype = gts.str();
        }
        cout << var << endl;
    }
    return 0;

}

