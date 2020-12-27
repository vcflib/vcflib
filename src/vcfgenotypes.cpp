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

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

    if (argc != 2) {
      cerr << "usage: " << argv[0] << " <vcf file>" << endl << endl
           << "report the genotypes for each sample, for each variant in the vcf file" << endl;
            cerr << R"(
Example:

      vcfgenotypes samples/sample.vcf

19      111     A       C       A,C     NA00001:A/A     NA00002:A/A     NA00003:A/C
19      112     A       G       A,G     NA00001:A/A     NA00002:A/A     NA00003:A/G
20      14370   G       A       G,A     NA00001:G/G     NA00002:G/A     NA00003:A/A
20      17330   T       A       T,A     NA00001:T/T     NA00002:T/A     NA00003:T/T
20      1110696 A       G,T     A,G,T   NA00001:G/T     NA00002:G/T     NA00003:T/T
20      1230237 T       .       T,.     NA00001:T/T     NA00002:T/T     NA00003:T/T
20      1234567 G       GA,GAC  G,GA,GAC        NA00001:G/GA    NA00002:G/GAC   NA00003:GA/GA
20      1235237 T       .       T,.     NA00001:T/T     NA00002:T/T     NA00003:./.
X       10      AC      A,ATG   AC,A,ATG        NA00001:AC      NA00002:AC/A    NA00003:AC/ATG
)";
      cerr << endl << "Type: statistics" << endl << endl;
        return 1;
    }

    string filename = argv[1];

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
    while (variantFile.getNextVariant(var)) {
        map<string, map<string, vector<string> > >::iterator s     = var.samples.begin();
        map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();

        cout << var.sequenceName << "\t"
             << var.position     << "\t"
             << var.ref          << "\t";
        var.printAlt(cout);     cout << "\t";
        var.printAlleles(cout); cout << "\t";

        for (; s != sEnd; ++s) {
            map<string, vector<string> >& sample = s->second;
            string& genotype = sample["GT"].front(); // XXX assumes we can only have one GT value
            map<int, int> gt = decomposeGenotype(genotype);
            // report the sample and it's genotype
            cout << s->first << ":";
            vector<string> x;
            for (map<int, int>::iterator g = gt.begin(); g != gt.end(); ++g) {
                for (int i = 0; i < g->second; ++i) {
                    if (g->first == NULL_ALLELE) {
                        x.push_back(".");
                    } else {
                        x.push_back(var.alleles[g->first]);
                    }
                }
            }
            cout << join(x, "/");
            cout << "\t";
        }
        cout << endl;
    }
    return 0;

}
