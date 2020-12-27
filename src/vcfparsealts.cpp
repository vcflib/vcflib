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
Alternate allele parsing method.  This method uses pairwise
alignment of REF and ALTs to determine component allelic primitives
for each alternate allele.

Usage: vcfparsealts <vcf file>

Example:

vcfparsealts samples/sample.vcf
##fileformat=VCFv4.0
(...)
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001 NA00002 NA00003
19      111     .       A       C       9.6     .       .       GT:HQ   0|0:10,10       0|0:10,10       0/1:3,3
 ( A :: 111 A -> A;  )  ( C :: 111 A -> C;  )
19      112     .       A       G       10      .       .       GT:HQ   0|0:10,10       0|0:10,10       0/1:3,3
 ( A :: 112 A -> A;  )  ( G :: 112 A -> G;  )
20      14370   rs6054257       G       A       29      PASS    AF=0.5;DP=14;NS=3;DB;H2 GT:GQ:DP:HQ     0|0:48:1:51,51     1|0:48:8:51,51  1/1:43:5:.,.
 ( A :: 14370 G -> A;  )  ( G :: 14370 G -> G;  )
20      17330   .       T       A       3       q10     AF=0.017;DP=11;NS=3     GT:GQ:DP:HQ     0|0:49:3:58,50     0|1:3:5:65,3    0/0:41:3:.,.
 ( A :: 17330 T -> A;  )  ( T :: 17330 T -> T;  )
20      1110696 rs6040355       A       G,T     67      PASS    AA=T;AF=0.333,0.667;DP=10;NS=2;DB       GT:GQ:DP:HQ        1|2:21:6:23,27  2|1:2:0:18,2    2/2:35:4:.,.
 ( A :: 1110696 A -> A;  )  ( G :: 1110696 A -> G;  )  ( T :: 1110696 A -> T;  )
20      1230237 .       T       .       47      PASS    AA=T;DP=13;NS=3 GT:GQ:DP:HQ     0|0:54:.:56,60  0|0:48:4:51,51     0/0:61:2:.,.
 ( . :: 1230237 T -> .;  )  ( T :: 1230237 T -> T;  )
20      1234567 microsat1       G       GA,GAC  50      PASS    AA=G;AC=3,1;AN=6;DP=9;NS=3      GT:GQ:DP  0/1:.:4  0/2:17:2        1/1:40:3
 ( G :: 1234567 G -> G;  )  ( GA :: 1234567 G -> G; 1234568  -> A;  )  ( GAC :: 1234567 G -> G; 1234568  -> AC;  )
20      1235237 .       T       .       0       .       .       GT      0/0     0|0     ./.
 ( . :: 1235237 T -> .;  )  ( T :: 1235237 T -> T;  )
X       10      rsTest  AC      A,ATG   10      PASS    .       GT      0       0/1     0|2
 ( A :: 10 A -> A; 11 C -> ;  )  ( AC :: 10 AC -> AC;  )  ( ATG :: 10 A -> A; 11  -> T; 11 C -> G;  )


Type: statistics

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
        map<string, vector<VariantAllele> > variants = var.parsedAlternates();
	cout << var << endl;
        for (map<string, vector<VariantAllele> >::iterator va = variants.begin(); va != variants.end(); ++va) {
            cout << " ( " << va->first << " :: ";
            vector<VariantAllele>& vars = va->second;
            vector<VariantAllele>::iterator g = vars.begin();
            for (; g != vars.end(); ++g) {
                cout << *g << "; ";
            }
            cout << " ) ";
        }
        cout << endl;
    }

    return 0;

}
