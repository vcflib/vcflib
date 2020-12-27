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
usage: vcfcat [file1] [file2] ... [fileN]

Concatenates VCF files

Type: transformation

      )";
    }
    exit(1);
  }

    if (argc == 1) {
        return 0;
    } else {
        for (int i = 1; i < argc; ++i) {
            VariantCallFile variantFile;
            string filename = argv[i];
            variantFile.open(filename);
            if (!variantFile.is_open()) {
                cerr << "could not open " << argv[i] << endl;
                return 1;
            }
            if (i == 1) {
                cout << variantFile.header << endl;
            }
            Variant var(variantFile);
            while (variantFile.getNextVariant(var)) {
                cout << var << endl;
            }
        }
    }

    return 0;

}
