#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcf;


int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "count the number of alternate alleles in heterozygous genotypes in all records in the vcf file" << endl
             << "outputs a count for each individual in the file" << endl;
        return 1;
    }

    string filename = argv[1];

    VariantCallFile variantFile;
    variantFile.open(filename);
    if (!variantFile.is_open()) {
        cerr << "could not open " << filename << endl;
        return 1;
    }

    unsigned int hetAlleleCount = 0;
    map<string, unsigned int> hetCounts;
    for (vector<string>::iterator s = variantFile.sampleNames.begin(); s != variantFile.sampleNames.end(); ++s) {
        hetCounts[*s] = 0;
    }

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        //cout << var << endl;
        for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
            string name = s->first;
            map<string, vector<string> >& sample = s->second;
            string& genotype = sample["GT"].front();
            vector<string> gt = split(genotype, "|/");
            int alt = 0;
            for (vector<string>::iterator g = gt.begin(); g != gt.end(); ++g) {
                if (*g != "0")
                    ++alt;
            }
            if (alt != gt.size()) {
                hetCounts[name] += alt;
                //hetAlleleCount += alt;
            }
        }
    }

    //cout << hetAlleleCount << endl;
    for (vector<string>::iterator s = variantFile.sampleNames.begin(); s != variantFile.sampleNames.end(); ++s) {
        cout << (s == variantFile.sampleNames.begin() ? "" : "\t") << *s;
    }
    cout << endl;
    for (vector<string>::iterator s = variantFile.sampleNames.begin(); s != variantFile.sampleNames.end(); ++s) {
        cout << (s == variantFile.sampleNames.begin() ? "" : "\t") << hetCounts[*s];
    }
    cout << endl;

    return 0;

}

