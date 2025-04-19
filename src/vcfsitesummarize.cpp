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
Summarize by site

Usage: vcfsitesummarize <vcf file>

Example:

vcfsitesummarize samples/sample.vcf

CHROM   POS     ID      REF     QUAL    FILTER  AA      AC      AF      AN      DP      NS      DB      H2
19      111     .       A       9.6     .                                                       0       0
19      112     .       A       10      .                                                       0       0
20      14370   rs6054257       G       29      PASS                    0.5             14      3       1 1
20      17330   .       T       3       q10                     0.017           11      3       0       0
20      1110696 rs6040355       A       67      PASS    T                               10      2       1 0
20      1230237 .       T       47      PASS    T                               13      3       0       0
20      1234567 microsat1       G       50      PASS    G                       6       9       3       0 0
20      1235237 .       T       0       .                                                       0       0
X       10      rsTest  AC      10      PASS


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

    // obtain all possible field names
    vector<string> infofields;
    vector<string> infoflags;

    for (map<string, VariantFieldType>::iterator i = variantFile.infoTypes.begin(); i != variantFile.infoTypes.end(); ++i) {
        if (variantFile.infoCounts[i->first] != ALLELE_NUMBER) {
            if (i->second == FIELD_BOOL) {
                infoflags.push_back(i->first);
            } else {
                infofields.push_back(i->first);
            }
        }
    }

    // write header

    // defaults
    cout << "CHROM\tPOS\tID\tREF\tQUAL\tFILTER";

    // configurable info field
    for (const auto& infoField : infofields) {
        cout << "\t" << infoField;
    }
    for (const auto& infoFlag : infoflags) {
        cout << "\t" << infoFlag;
    }
    cout << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

	cout << var.sequenceName << "\t"
	     << var.position << "\t"
	     << var.id << "\t"
	     << var.ref << "\t"
	     << var.quality << "\t"
	     << var.filter;

	for (const auto& name : infofields) {
	    vector<string> value;
	    map<string, vector<string> >::iterator f = var.info.find(name);
	    if (f != var.info.end()) {
            value = f->second;
            if (value.size() == 1) {
                cout << "\t" << value.front();
            } else {
                cout << "\t"; // null
            }
	    } else {
            cout << "\t"; // null
	    }
	}

	for (const auto& name : infoflags) {
	    string value;
	    map<string, bool>::iterator f = var.infoFlags.find(name);
	    cout << "\t";
	    if (f != var.infoFlags.end()) {
            cout << 1;
	    } else {
            cout << 0;
	    }
	}

	cout << endl;

    }

    return 0;

}
