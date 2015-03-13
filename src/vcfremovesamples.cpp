#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcflib;

// remove elements in a from b
template<class T>
vector<T> removeElems(vector<T>& a, vector<T>& b) {
    map<T, bool> inA;
    map<T, bool> inAB;
    for (typename vector<T>::iterator i = a.begin(); i != a.end(); ++i) {
        inA[*i] = true;
    }
    for (typename vector<T>::iterator i = b.begin(); i != b.end(); ++i) {
        if (inA.find(*i) == inA.end()) {
            inAB[*i] = true;
        }
    }
    vector<T> aNb;
    for (typename map<T, bool>::iterator i = inAB.begin(); i != inAB.end(); ++i) {
        aNb.push_back(i->first);
    }
    return aNb;
}

int main(int argc, char** argv) {

    if (argc < 3) {
        cerr << "usage: " << argv[0] << " <vcf file> [SAMPLE1] [SAMPLE2] ..." << endl
             << "outputs each record in the vcf file, removing samples listed on the command line" << endl;
        return 1;
    }

    string filename = argv[1];

    vector<string> samplesToRemove;
    for (int i = 2; i < argc; ++i) {
        samplesToRemove.push_back(argv[i]);
    }

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

    vector<string> samplesToKeep = removeElems(samplesToRemove, variantFile.sampleNames);

    // update sample list in header
    variantFile.updateSamples(samplesToKeep);

    // and restrict the output sample names in the variant to those we are keeping
    var.setOutputSampleNames(samplesToKeep);
    
    // write the new header
    cout << variantFile.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        cout << var << endl;
    }

    return 0;

}

