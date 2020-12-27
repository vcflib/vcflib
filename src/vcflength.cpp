/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include "convert.h"
#include <vector>

using namespace std;
using namespace vcflib;

int main(int argc, char** argv) {
if (argc == 2) {
  string h_flag = argv[1];
  if (h_flag == "-h" || h_flag == "--help") {
    cerr << R"(
Add length info field

Usage: vcflength

Example:

vcflength samples/sample.vcf
##fileformat=VCFv4.0
(...)
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001 NA00002 NA00003
19      111     .       A       C       9.6     .       length=0;length.alt=1;length.ref=1      GT:HQ   0|0:10,10  0|0:10,10       0/1:3,3
19      112     .       A       G       10      .       length=0;length.alt=1;length.ref=1      GT:HQ   0|0:10,10  0|0:10,10       0/1:3,3
20      14370   rs6054257       G       A       29      PASS    AF=0.5;DP=14;NS=3;length=0;length.alt=1;length.ref=1;DB;H2 GT:GQ:DP:HQ     0|0:48:1:51,51  1|0:48:8:51,51  1/1:43:5:.,.
20      17330   .       T       A       3       q10     AF=0.017;DP=11;NS=3;length=0;length.alt=1;length.ref=1     GT:GQ:DP:HQ     0|0:49:3:58,50  0|1:3:5:65,3    0/0:41:3:.,.
20      1110696 rs6040355       A       G,T     67      PASS    AA=T;AF=0.333,0.667;DP=10;NS=2;length=0,0;length.alt=1,1;length.ref=1;DB   GT:GQ:DP:HQ     1|2:21:6:23,27  2|1:2:0:18,2    2/2:35:4:.,.
20      1230237 .       T       .       47      PASS    AA=T;DP=13;NS=3;length=0;length.alt=1;length.ref=1GT:GQ:DP:HQ      0|0:54:.:56,60  0|0:48:4:51,51  0/0:61:2:.,.
20      1234567 microsat1       G       GA,GAC  50      PASS    AA=G;AC=3,1;AN=6;DP=9;NS=3;length=1,2;length.alt=2,3;length.ref=1  GT:GQ:DP        0/1:.:4 0/2:17:2        1/1:40:3
20      1235237 .       T       .       0       .       length=0;length.alt=1;length.ref=1      GT      0/00|0     ./.
X       10      rsTest  AC      A,ATG   10      PASS    length=-1,1;length.alt=1,3;length.ref=2 GT      0 0/1      0|2

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

    variantFile.addHeaderLine("##INFO=<ID=length,Number=A,Type=Integer,Description=\"length(ALT) - length(REF) for each ALT\">");
    variantFile.addHeaderLine("##INFO=<ID=length.ref,Number=1,Type=Integer,Description=\"length(REF)\">");
    variantFile.addHeaderLine("##INFO=<ID=length.alt,Number=A,Type=Integer,Description=\"length(ALT) for each ALT\">");
    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        vector<string>& lengths = var.info["length"];
        lengths.clear();
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            lengths.push_back(convert((int) a->size() - (int) var.ref.size()));
        }
        vector<string>& lengthsRef = var.info["length.ref"];
        lengthsRef.clear();
        lengthsRef.push_back(convert(var.ref.size()));
        vector<string>& lengthsAlt = var.info["length.alt"];
        lengthsAlt.clear();
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            lengthsAlt.push_back(convert((int) a->size()));
        }
        cout << var << endl;
    }

    return 0;

}
