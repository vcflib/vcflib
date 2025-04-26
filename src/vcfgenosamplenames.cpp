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
Get samplenames

Usage: vcfgenosamplenames

Example:

vcfsamplenames samples/sample.vcf

NA00001
NA00002
NA00003


Type: transformation

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

    variantFile.addHeaderLine("##FORMAT=<ID=SN,Number=1,Type=String,Description=\"The name of the sample.\">");

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        var.format.push_back("SN");
        for (auto& s : var.samples) {
            s.second["SN"].clear();
            s.second["SN"].push_back(s.first);
        }
        cout << var << endl;
    }

    return 0;

}
