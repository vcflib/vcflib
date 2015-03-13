#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "outputs a VCF stream in which 'long' non-complex"
             << "alleles have their position corrected." << endl
             << "assumes that VCF records can't overlap 5'->3'" << endl;
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

    Variant var(variantFile);

    // write the new header
    cout << variantFile.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        // if we just have one parsed alternate (non-complex case)
        map<string, vector<VariantAllele> > parsedAlts = var.parsedAlternates(true, true); // use mnps, and previous for indels
        // but the alt string is long
        //cerr << var.alt.size() << " " << parsedAlts.size() << endl;
        if (var.alt.size() == 1 && parsedAlts.size() > 1) {
            string& alternate = var.alt.front();
            vector<VariantAllele>& vs = parsedAlts[alternate];
            vector<VariantAllele> valleles;
            for (vector<VariantAllele>::iterator a = vs.begin(); a != vs.end(); ++a) {
                if (a->ref != a->alt) {
                    valleles.push_back(*a);
                }
            }
            if (valleles.size() == 1) {
                // do we have extra sequence hanging around?
                VariantAllele& varallele = valleles.front();
                if (vs.front().ref == vs.front().alt) {
                    var.position = varallele.position;
                    var.ref = var.ref.substr(vs.front().ref.size(), varallele.ref.size());
                    var.alt.front() = varallele.alt;
                }
            }
        }
        cout << var << endl;
    }

    return 0;

}

