#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>
#include <set>

using namespace std;
using namespace vcflib;

// adds non-overlapping info fields from varB to varA
void addInfo(Variant& varA, Variant& varB) {
    for (map<string, vector<string> >::iterator i = varB.info.begin(); i != varB.info.end(); ++i) {
        if (varA.info.find(i->first) == varA.info.end()) {
            varA.info[i->first] = i->second;
        }
    }
}

int main(int argc, char** argv) {

    if (argc != 3) {
        cerr << "usage: " << argv[0] << " <vcf file> <vcf file>" << endl
             << "Adds info fields from the second file which are not present in the first vcf file." << endl;
        return 1;
    }

    string filenameA = argv[1];
    string filenameB = argv[2];

    if (filenameA == filenameB) {
        cerr << "it won't help to add info data from the same file!" << endl;
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

    // while the first file doesn't match the second positionally,
    // step forward, annotating each genotype record with an empty genotype
    // when the two match, iterate through the genotypes from the first file
    // and get the genotypes reported in the second file
    
    variantFileA.getNextVariant(varA);
    variantFileB.getNextVariant(varB);
    
    variantFileA.header = unionInfoHeaderLines(variantFileA.header, variantFileB.header);
    
    cout << variantFileA.header << endl;

    do {

        while (!variantFileB.done()
               && (varB.sequenceName < varA.sequenceName
                   || (varB.sequenceName == varA.sequenceName && varB.position < varA.position))
            ) {
            variantFileB.getNextVariant(varB);
        }

        while (!variantFileA.done()
               && (varA.sequenceName < varB.sequenceName
                   || (varA.sequenceName == varB.sequenceName && varA.position < varB.position))
            ) {
            cout << varA << endl;
            variantFileA.getNextVariant(varA);
        }

        while (!variantFileB.done()
               && (varB.sequenceName < varA.sequenceName
                   || (varB.sequenceName == varA.sequenceName && varB.position < varA.position))
            ) {
            variantFileB.getNextVariant(varB);
        }

        while (!variantFileA.done() && varA.sequenceName == varB.sequenceName && varA.position == varB.position) {
            addInfo(varA, varB);
            cout << varA << endl;
            variantFileA.getNextVariant(varA);
            variantFileB.getNextVariant(varB);
        }
        
    } while (!variantFileA.done() && !variantFileB.done());

    if (!variantFileA.done()) {
        cout << varA << endl;
        while (variantFileA.getNextVariant(varA)) {
            cout << varA << endl;
        }
    }

    return 0;

}

