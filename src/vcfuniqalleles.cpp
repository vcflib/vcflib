/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include <set>

using namespace std;
using namespace vcflib;

int main(int argc, char** argv) {
if (argc == 2) {
  string h_flag = argv[1];
  if (h_flag == "-h" || h_flag == "--help") {
      cerr << R"(

List unique alleles For each record, remove any duplicate alternate
alleles that may have resulted from merging separate VCF files.

Usage: vcfuniqalleles <vcf file>

Type: filter

      )";
      exit(1);
    }
  }

    VariantCallFile variantFile;

    if (argc > 1) {
        string filename = argv[1];
        variantFile.open(filename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        set<string> alleles;
        vector<string> alleles_to_remove;
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            if (*a != var.ref) {
                if (alleles.find(*a) == alleles.end()) {
                    alleles.insert(*a);
                } else {
                    alleles_to_remove.push_back(*a);
                }
            } else {
                alleles_to_remove.push_back(*a); // same as ref
            }
        }
        for (vector<string>::iterator a = alleles_to_remove.begin(); a != alleles_to_remove.end(); ++a) {
            cerr << "removing " << *a << " from " << var.sequenceName << ":" << var.position << endl;
            var.removeAlt(*a);
        }
        cout << var << endl;
    }

    return 0;

}
