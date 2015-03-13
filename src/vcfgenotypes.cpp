#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "report the genotypes for each sample, for each variant in the vcf file" << endl;
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

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        map<string, map<string, vector<string> > >::iterator s     = var.samples.begin(); 
        map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();
        
        cout << var.sequenceName << "\t"
             << var.position     << "\t"
             << var.ref          << "\t";
        var.printAlt(cout);     cout << "\t"; 
        var.printAlleles(cout); cout << "\t"; 
        
        for (; s != sEnd; ++s) {
            map<string, vector<string> >& sample = s->second;
            string& genotype = sample["GT"].front(); // XXX assumes we can only have one GT value
            vector<string> gt = split(genotype, "|/");
            
            // report the sample and it's genotype
            cout << s->first << ":";
            for (vector<string>::iterator g = gt.begin(); g != gt.end(); ++g) {
                if (g->c_str() == ".") {
                    cout << ".";
                } else {
                    int index = atoi(g->c_str());
                    cout << var.alleles[index];
                }
                if (g != (gt.end()-1)) cout << "/";
            }
            cout << "\t";
        }
        cout << endl;
    }
    return 0;

}

