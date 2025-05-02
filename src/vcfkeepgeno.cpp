/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"

#include <string>
#include <iostream>
#include <set>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

  if (argc == 2) {
    string h_flag = argv[1];
    if (h_flag == "-h" || h_flag == "--help") {
      cerr << "usage: " << argv[0] << " <vcf file> [FIELD1] [FIELD2] ..." << endl << endl
           << "Reduce file size by removing FORMAT fields not listed "
           << "on the command line from sample specifications in the output"
           << endl;
      cerr << endl << "Type: transformation" << endl << endl;
      return 1;
    }
  }

    string filename = argv[1];

    vector<string> newFormat;
    set<string> fieldsToKeep;
    for (int i = 2; i < argc; ++i) {
        fieldsToKeep.insert(argv[i]);
        newFormat.push_back(argv[i]);
    }

    VariantCallFile variantFile;
    if (filename == "-") {
        variantFile.open(std::cin);
    } else {
        variantFile.open(filename);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    Variant var(variantFile);

    vector<string> formatIds = variantFile.formatIds();
    for (const auto& i : formatIds) {
        if (!fieldsToKeep.count(i)) {
            variantFile.removeGenoHeaderLine(i);
        }
    }

    // write the header
    cout << variantFile.header << endl;

    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        var.format = newFormat;
        cout << var << endl;
    }

    return 0;

}
