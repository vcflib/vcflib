#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "outputs the het/hom ratio for each individual in the file" << endl;
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
        cerr << "could not open " << filename << endl;
        return 1;
    }

    map<string, unsigned int> hetCounts;
    map<string, unsigned int> homCounts;
    for (vector<string>::iterator s = variantFile.sampleNames.begin(); s != variantFile.sampleNames.end(); ++s) {
        hetCounts[*s] = 0;
        homCounts[*s] = 0;
    }

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        //cout << var << endl;
        for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
            string name = s->first;
            map<string, vector<string> >& sample = s->second;
            string& gt = sample["GT"].front();
            map<int, int> genotype = decomposeGenotype(gt);
            if (isHet(genotype)) {
                ++hetCounts[name];
            } else if (isHomNonRef(genotype)) {
                ++homCounts[name];
            }
        }
    }

    for (vector<string>::iterator s = variantFile.sampleNames.begin(); s != variantFile.sampleNames.end(); ++s) {
        cout << (s == variantFile.sampleNames.begin() ? "" : "\t") << *s;
    }
    cout << endl;
    for (vector<string>::iterator s = variantFile.sampleNames.begin(); s != variantFile.sampleNames.end(); ++s) {
        cout << (s == variantFile.sampleNames.begin() ? "" : "\t") << (double) hetCounts[*s] / (double) homCounts[*s];
    }
    cout << endl;

    return 0;

}

