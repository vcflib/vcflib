#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
using namespace vcf;

void stripAberrant(Variant& var) {
    map<string, map<string, string> >::iterator s = var.samples.begin();
    while (s != var.samples.end()) {
        map<string, string>& sample = s->second;
        map<string, int> genotype = decomposeGenotype(sample["GT"]);
        int altobs = 0;
        int refobs = 0;
        convert(sample["RA"], refobs);
        convert(sample["AA"], altobs);
        if (isHomRef(genotype) && altobs > 0
            || isHomNonRef(genotype) && refobs > 0) {
            var.samples.erase(s++);
        } else {
            ++s;
        }
    }
}

int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "strips samples which are homozygous but have observations implying heterozygosity" << endl;
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

    // TODO check if AC is present
    // ensure that AC is listed as an info field
    string line = "##filter=\"removed homozygous genotypes which have observations implying heterozygosity\">";
    variantFile.addHeaderLine(line);

    // write the new header
    cout << variantFile.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        stripAberrant(var);
        cout << var << endl;
    }

    return 0;

}

