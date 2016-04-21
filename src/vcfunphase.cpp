#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "unphases and sorts the genotypes in the file" << endl;
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
            vector<string> gt = split(genotype, "|/");
            // now let's sort the genotype
            vector<int> gti;
            for (vector<string>::iterator g = gt.begin(); g != gt.end(); ++g) {
                if (*g == ".") {
                    gti.push_back(-1);
                } else {
                    gti.push_back(atoi(g->c_str()));
                }
            }
            std::sort(gti.begin(), gti.end());
            stringstream gts;
            for (vector<int>::iterator g = gti.begin(); g != gti.end(); ++g) {
                if (g != gti.begin()) {
                    gts << "/";
                }
                if (*g == -1) {
                    gts << ".";
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

