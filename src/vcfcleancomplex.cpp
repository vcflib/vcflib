/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020-2024      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "legacy.h"

#include <string>
#include <sstream>
#include <iostream>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

  if (argc == 2) {
    string h_flag = argv[1];
    if (h_flag == "-h" || h_flag == "--help") {
      cerr << "usage: " << argv[0] << " <vcf file>" << endl << endl
           << "Removes reference-matching sequence from complex alleles and adjusts records to reflect positional change." << endl << endl
           << "Generate a VCF stream in which 'long' non-complex"
           << "alleles have their position corrected." << endl
           << "assumes that VCF records can't overlap 5'->3'" << endl;
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

    VariantLegacy var(variantFile);

    // write the new header
    cout << variantFile.header << endl;

    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        // if we just have one parsed alternate (non-complex case)
        map<string, vector<VariantAllele> > parsedAlts = var.legacy_parsedAlternates(true, true); // use mnps, and previous for indels
        // but the alt string is long
        //cerr << var.alt.size() << " " << parsedAlts.size() << endl;
        if (var.alt.size() == 1 && parsedAlts.size() > 1) {
            string& alternate = var.alt.front();
            vector<VariantAllele>& vs = parsedAlts[alternate];
            vector<VariantAllele> valleles;
            for (const auto& a : vs) {
                if (a.ref != a.alt) {
                    valleles.push_back(a);
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
