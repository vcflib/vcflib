/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
using namespace vcflib;

bool isTransition(const string& ref, const string& alt) {
    if (((ref == "A" && alt == "G") || (ref == "G" && alt == "A")) ||
        ((ref == "C" && alt == "T") || (ref == "T" && alt == "C"))) {
        return true;
    } else {
        return false;
    }
}

bool hasTransition(const Variant& var) {
    for (const auto& alt : var.alt) {
        if (isTransition(var.ref, alt)) {
            return true;
        }
    }
    return false;
}

bool hasTransversion(const Variant& var) {
    for (const auto& alt : var.alt) {
        if (!isTransition(var.ref, alt)) {
            return true;
        }
    }
    return false;
}

bool hasInsertion(Variant& var) {
    string& ref = var.ref;
    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
        string& alt = *a;
        if (ref.size() < alt.size()) {
            return true;
        }
    }
    return false;
}

bool hasDeletion(Variant& var) {
    for (const auto& alt : var.alt) {
        if (var.ref.size() > alt.size()) {
            return true;
        }
    }
    return false;
}

bool hasMNP(Variant& var) {
    for (const auto& alt : var.alt) {
        if (var.ref.size() > 1 && alt.size() == var.ref.size()) {
            return true;
        }
    }
    return false;
}

bool hasSNP(Variant& var) {
    for (const auto& alt : var.alt) {
        if (var.ref.size() == 1 && alt.size() == 1) {
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv) {

  if (argc == 2) {
    string h_flag = argv[1];
    if (h_flag == "-h" || h_flag == "--help") {
      cerr << R"(
usage: vcfclassify <vcf file>

Creates a new VCF where each variant is tagged by allele class: snp,
ts/tv, indel, mnp

Type: transformation

      )";
      exit(1);
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

    string line;
    line = "##INFO=<ID=SNP,Number=0,Type=Flag,Description=\"SNP allele\">";
    variantFile.addHeaderLine(line);
    line = "##INFO=<ID=TS,Number=0,Type=Flag,Description=\"transition SNP\">";
    variantFile.addHeaderLine(line);
    line = "##INFO=<ID=TV,Number=0,Type=Flag,Description=\"transversion SNP\">";
    variantFile.addHeaderLine(line);
    line = "##INFO=<ID=INS,Number=0,Type=Flag,Description=\"insertion allele\">";
    variantFile.addHeaderLine(line);
    line = "##INFO=<ID=DEL,Number=0,Type=Flag,Description=\"deletion allele\">";
    variantFile.addHeaderLine(line);
    line = "##INFO=<ID=MNP,Number=0,Type=Flag,Description=\"MNP allele\">";
    variantFile.addHeaderLine(line);
    // TODO handle lengths at poly-allelic sites
    //line = "##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"allele length\">";
    //variantFile.addHeaderLine(line);

    // write the new header
    cout << variantFile.header << endl;


    while (variantFile.getNextVariant(var)) {

        if (hasSNP(var)) {
            var.infoFlags["SNP"] = true;
        }

        if (hasTransition(var)) {
            var.infoFlags["TS"] = true;
        }

        if (hasTransversion(var)) {
            var.infoFlags["TV"] = true;
        }

        if (hasInsertion(var)) {
            var.infoFlags["INS"] = true;
        }

        if (hasDeletion(var)) {
            var.infoFlags["DEL"] = true;
        }

        if (hasMNP(var)) {
            var.infoFlags["MNP"] = true;
        }

        cout << var << endl;
    }

    return 0;

}
