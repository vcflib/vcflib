/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "convert.h"

#include "Variant.h"
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

using namespace std;
using namespace vcflib;

void stripAberrant(Variant& var) {
    map<string, map<string, vector<string> > >::iterator s = var.samples.begin();
    while (s != var.samples.end()) {
        map<string, vector<string> >& sample = s->second;
        map<int, int> genotype = decomposeGenotype(sample["GT"].front());
        int refobs = 0;
        convert(sample["RO"].front(), refobs);
        bool eraseSample = false;
        if (isHomNonRef(genotype) && refobs > 0) {
            eraseSample = true;
        } else if (isHomRef(genotype)) {
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                int alleleIndex = var.altAlleleIndexes[*a];
                int altobs = 0;
                convert(sample["AO"].at(alleleIndex), altobs);
                if (altobs > 0) {
                    eraseSample = true;
                    break;
                }
            }
        }
        if (eraseSample) {
            s = var.samples.erase(s);
        } else {
            ++s;
        }
    }
}

int main(int argc, char** argv) {

  if (argc == 2) {
    string h_flag = argv[1];
    if (h_flag == "-h" || h_flag == "--help") {
      cerr << "usage: " << argv[0] << " <vcf file>" << endl << endl
           << "strips samples which are homozygous but have observations implying heterozygosity." << endl
           << "Remove samples for which the reported genotype (GT) and observation counts disagree (AO, RO)." << endl;
      cerr << endl << "Type: transformation" << endl << endl;
      return 1;
    }
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
