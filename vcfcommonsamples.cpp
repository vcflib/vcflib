#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcf;

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

    map<string, bool> samplesInB;
    // get the samples from the second file
    for (vector<string>::iterator s = variantFileB.sampleNames.begin();
            s != variantFileB.sampleNames.end(); ++s) {
        samplesInB[*s] = true;
    }

    // update sample list in header
    variantFileA.updateSamples(variantFileB.sampleNames);

    // and restrict the output sample names in the variant to those we are keeping
    varA.setOutputSampleNames(variantFileB.sampleNames);
    
    // write the new header
    cout << variantFileA.header << endl;
    
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFileA.getNextVariant(varA)) {
        cout << varA << endl;
    }

    return 0;

}

