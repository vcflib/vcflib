#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcflib;

template<class T>
vector<T> intersection(vector<T>& a, vector<T>& b) {
    map<T, bool> inA;
    map<T, bool> inAB;
    for (typename vector<T>::iterator i = a.begin(); i != a.end(); ++i) {
        inA[*i] = true;
    }
    for (typename vector<T>::iterator i = b.begin(); i != b.end(); ++i) {
        if (inA.find(*i) != inA.end()) {
            inAB[*i] = true;
        }
    }
    vector<T> aIb;
    for (typename map<T, bool>::iterator i = inAB.begin(); i != inAB.end(); ++i) {
        aIb.push_back(i->first);
    }
    return aIb;
}

int main(int argc, char** argv) {

    if (argc != 3) {
        cerr << "usage: " << argv[0] << " <vcf file> <vcf file>" << endl
             << "outputs each record in the first file, removing samples not present in the second" << endl;
        return 1;
    }

    string filenameA = argv[1];
    string filenameB = argv[2];

    if (filenameA == filenameB) {
        cerr << "you're just spinning your wheels matching the samples in "
            << filenameA << " to the samples in " << filenameB << endl;
        return 1;
    }

    VariantCallFile variantFileA;
    if (filenameA == "-") {
        variantFileA.open(std::cin);
    } else {
        variantFileA.open(filenameA);
    }

    VariantCallFile variantFileB;
    if (filenameB == "-") {
        variantFileB.open(std::cin);
    } else {
        variantFileB.open(filenameB);
    }

    if (!variantFileA.is_open() || !variantFileB.is_open()) {
        return 1;
    }

    Variant varA(variantFileA);
    Variant varB(variantFileB);

    vector<string> commonSamples = intersection(variantFileA.sampleNames, variantFileB.sampleNames);

    // update sample list in header
    variantFileA.updateSamples(commonSamples);

    // and restrict the output sample names in the variant to those we are keeping
    varA.setOutputSampleNames(commonSamples);
 
    // write the new header
    cout << variantFileA.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFileA.getNextVariant(varA)) {
        cout << varA << endl;
    }

    return 0;

}

