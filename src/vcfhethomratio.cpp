/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"

#include <string>
#include <iostream>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

  if (argc == 2) {
    string h_flag = argv[1];
    if (h_flag == "-h" || h_flag == "--help") {
      cerr << "usage: " << argv[0] << " <vcf file>" << endl << endl
           << "Generates the het/hom ratio for each individual in the file" << endl;
      cerr << endl << "Type: metrics" << endl << endl;
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

    map<string, unsigned int> hetCounts;
    map<string, unsigned int> homCounts;
    for (const auto& s : variantFile.sampleNames) {
        hetCounts[s] = 0;
        homCounts[s] = 0;
    }

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        //cout << var << endl;
        for (auto& [name, sample] : var.samples) {
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
