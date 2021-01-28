/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"

using namespace std;
using namespace vcflib;

int main(int argc, char** argv) {

  if (argc == 2) {
    string h_flag = argv[1];
    if (h_flag == "-h" || h_flag == "--help") {
      cerr << R"(

Usage: vcfuniq <vcf file>

List unique genotypes.  Like GNU uniq, but for VCF records.  Remove
records which have the same position, ref, and alt as the previous
record.

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

    string lastsn;
    long int lastpos;
    string lastref;
    vector<string> lastalt;

    variantFile.parseSamples = false;
    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        if (!lastsn.empty()
            && (lastsn == var.sequenceName
                && lastpos == var.position
                && lastref == var.ref
                && lastalt == var.alt)) {
            continue;
        } else {
            lastsn = var.sequenceName;
            lastpos = var.position;
            lastref = var.ref;
            lastalt = var.alt;
            cout << var.originalLine << endl;
        }
    }

    return 0;

}
