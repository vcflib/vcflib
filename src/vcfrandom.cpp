/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include <sstream>
#include <stdlib.h>
#include <time.h>
#include "Variant.h"
#include <cmath>

using namespace std;
using namespace vcflib;

int main(int argc, char** argv) {

  if (argc == 2) {
    string h_flag = argv[1];

    if (argc == 2 && (h_flag == "-h" || h_flag == "--help")) {
      cerr << R"(
Generate a random VCF file

Usage: vcfrandom

Example:

    vcfrandom

##fileformat=VCFv4.0
##source=vcfrandom
##reference=/d2/data/references/build_37/human_reference_v37.fa
##phasing=none
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  bill
one     1       .       G       G,A     100     .       DP=83   GT:DP   0/1:1
one     2       .       G       G,A     100     .       DP=3    GT:DP   0/1:49
one     3       .       G       C,T     100     .       DP=5    GT:DP   0/1:12
one     4       .       C       G,T     100     .       DP=51   GT:DP   0/1:60
one     5       .       A       T,A     100     .       DP=31   GT:DP   0/1:89
one     6       .       T       T,A     100     .       DP=56   GT:DP   0/1:60
one     7       .       T       A,C     100     .       DP=78   GT:DP   0/1:75
one     8       .       T       G,A     100     .       DP=73   GT:DP   0/1:78
one     9       .       C       C,G     100     .       DP=42   GT:DP   0/1:67


Type: statistics

      )";
      exit(1);
    }
  }


    VariantCallFile variantFile;

    stringstream headerss;
    headerss << "##fileformat=VCFv4.0" << endl
             << "##source=vcfrandom" << endl
             << "##reference=/d2/data/references/build_37/human_reference_v37.fa" << endl
             << "##phasing=none" << endl
             << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << endl
             << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">" << endl
             << "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">" << endl
             << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl
             << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">" << endl
             << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
             << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype\">" << endl
             << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl
             << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbill";

    string header = headerss.str();
    variantFile.openForOutput(header);

    cout << variantFile.header << endl;

    srand(time(NULL));

    vector<string> atgc;
    atgc.push_back("A");
    atgc.push_back("T");
    atgc.push_back("G");
    atgc.push_back("C");

    for (int i = 1; i < 10; ++i) {
        Variant var(variantFile);
        var.sequenceName = "one";
        var.id = ".";
        var.filter = ".";
        var.ref = atgc.at(rand() % 4);
        var.quality = 100;
        stringstream s;
        s << rand() % 100;
        var.info["DP"].push_back(s.str());
        var.format.push_back("GT");
        var.format.push_back("DP");
        var.position = i;
        for (vector<string>::iterator s = var.sampleNames.begin(); s != var.sampleNames.end(); ++s) {
            string& name = *s;
            var.alt.clear();
            var.alt.push_back(atgc.at(rand() % 4));
            var.alt.push_back(atgc.at(rand() % 4));
            var.samples[name]["GT"].push_back("0/1");
            stringstream dp;
            dp << floor(rand() % 100);
            var.samples[name]["DP"].push_back(dp.str());
        }
        cout << var << endl;
    }

    return 0;

}
