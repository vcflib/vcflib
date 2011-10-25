#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
using namespace vcf;


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
    cout << variantFile.header;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        // if we just have one parsed alternate (non-complex case)
        map<string, vector<VariantAllele> > parsedAlts = var.parsedAlternates();
        // but the alt string is long
        if (var.alt.size() == 1 && parsedAlts.size() == 2) {
            string& alternate = var.alt.front();
            if (parsedAlts[alternate].size() == 1) {
                // do we have extra sequence hanging around?
                VariantAllele& varallele = parsedAlts[alternate].front();
                // for deletions and insertions, we have to keep a leading base
                // but the variant allele doesn't have these
                if (varallele.ref.size() == varallele.alt.size()) {
                    if (varallele.position != var.position) {
                        var.ref = varallele.ref;
                        var.alt.front() = varallele.alt;
                        var.position = varallele.position;
                    }
                } else if (varallele.ref.size() < varallele.alt.size()) {
                    if (varallele.position != var.position + 1) {
                        // TODO unhandled
                    }
                } else if (varallele.ref.size() < varallele.alt.size()) {
                    if (varallele.position != var.position) {
                        // TODO unhandled
                    }
                }
            }
        }
        cout << var << endl;
    }

    return 0;

}

