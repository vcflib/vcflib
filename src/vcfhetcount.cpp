/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>
#include <getopt.h>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

  if (argc == 2) {
    string h_flag = argv[1];
    if (h_flag == "-h" || h_flag == "--help") {
      cerr << "usage: " << argv[0] << " <vcf file>" << endl << endl
           << "Calculate the heterozygosity rate: " << endl
           << "count the number of alternate alleles in heterozygous genotypes in all records in the vcf file" << endl
           << "outputs a count for each individual in the file" << endl;
      cerr << endl << "Type: metrics" << endl << endl;
        return 1;
    }
  }


    string inputFilename;
    VariantCallFile variantFile;
    if (optind == argc - 1) {
        inputFilename = argv[optind];
        variantFile.open(inputFilename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    unsigned int hetAlleleCount = 0;
    map<string, unsigned int> hetCounts;
    for (const auto& sampleName : variantFile.sampleNames) {
        hetCounts[sampleName] = 0;
    }

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        //cout << var << endl;
        for (auto& [name, sample] : var.samples) {
            string& genotype = sample["GT"].front();
            vector<string> gt = split(genotype, "|/");
            int alt = 0;
            for (const auto& g : gt) {
                if (g != "0")
                    ++alt;
            }
            if (alt != gt.size()) {
                hetCounts[name] += alt;
                //hetAlleleCount += alt;
            }
        }
    }

    //cout << hetAlleleCount << endl;
    for (vector<string>::iterator s = variantFile.sampleNames.begin(); s != variantFile.sampleNames.end(); ++s) {
        cout << (s == variantFile.sampleNames.begin() ? "" : "\t") << *s;
    }
    cout << endl;
    for (vector<string>::iterator s = variantFile.sampleNames.begin(); s != variantFile.sampleNames.end(); ++s) {
        cout << (s == variantFile.sampleNames.begin() ? "" : "\t") << hetCounts[*s];
    }
    cout << endl;

    return 0;

}
