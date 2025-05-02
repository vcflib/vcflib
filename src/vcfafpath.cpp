/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "convert.h"

#include "Variant.h"
#include <vector>
#include <map>

#include "join.h"

using namespace std;
using namespace vcflib;

int main(int argc, char** argv) {
  if (argc == 2) {
    string h_flag = argv[1];

    if (argc == 2 && (h_flag == "-h" || h_flag == "--help")) {
      cerr << R"(
Display genotype paths

Usage: vcfafpath <vcf file>

Example:

    vcfafpath samples/scaffold612.vcf

```

T -> A
A -> G
T -> C
C -> A
C -> T
A -> G
T -> C
G -> C
C -> CAGA
A -> G
```


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

    //cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        //cout << var << endl;
        double afref = 1;
        map<double, vector<string> > allelesByAf;
        vector<double> afd;
        vector<string>& afstr = var.info["AF"];
        for (vector<string>::iterator af = afstr.begin(); af != afstr.end(); ++af) {
            double r; convert(*af, r);
            afd.push_back(r);
        }
        vector<double>::iterator af = afd.begin();
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a, ++af) {
            afref -= *af;
            allelesByAf[*af].push_back(*a);
        }
        cout << var.ref;
        for (map<double, vector<string> >::reverse_iterator a = allelesByAf.rbegin(); a != allelesByAf.rend(); ++a) {
            cout << " -> " << join(a->second, ", ");
        }
        cout << endl;
    }

    return 0;

}
