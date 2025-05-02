/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2022 Erik Garrison
    Copyright © 2020-2022 Pjotr Prins
    Copyright © 2015-2022 Travis Collier

    This software is published under the MIT License. See the LICENSE file.
*/


#include "Variant.h"

#include <getopt.h>
#include <string>
#include <iostream>

using namespace std;
using namespace vcflib;

void printSummary(char** argv) {
    cout << "usage: " << argv[0] << " [options] <vcf file>";
    std::string text = R"(

Makes the FORMAT for each variant line the same (uses all the FORMAT
fields described in the header). Fills out per-sample fields to match
FORMAT. Expands GT values of '.' with number of alleles based on
ploidy (eg: './.' for dipolid).

options:

    -p, --ploidy N   the polidy of missing/null GT fields (default=2)
    -L, --expand_GL  fill in missing GL fields with 0 values (eg: 0,0,0 for diploid 2 alleles)

Type: transformation
)";
    cout << text;
    exit(0);
}

int nChoosek(int n, int k){
    if (k > n) return 0;
    if (k*2 > n) k = n-k;
    if (k == 0) return 1;
    int result = n;
    for( int i = 2; i <= k; ++i ){
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

int main(int argc, char** argv) {
    int ploidy = 2;
    bool expand_GL = false;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"expand_GL", no_argument, 0, 'L'},
                {"ploidy", required_argument, 0, 'p'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        c = getopt_long (argc, argv, "hp:L", long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {

        case 'p':
            ploidy = atoi(optarg);
            break;
        case 'L':
            expand_GL = true;
            break;
        case 'h':
        default:
            printSummary(argv);
            break;
        }
    }

    VariantCallFile variantFile;

    variantFile.open(std::cin);

    if (!variantFile.is_open()) {
        return 1;
    }

    cout << variantFile.header << endl;

    vector<string> null_val_vector;
    null_val_vector.push_back(".");

    vector<string> null_GTval_vector;
    null_GTval_vector.push_back(".");
    for( int i=1; i<ploidy; ++i ){
        null_GTval_vector.back() += "/.";
    }

    vector<string> null_GLval_vector;
    null_GLval_vector.push_back(".");

    Variant var(variantFile);

    vector<string> out_formats;
    out_formats.push_back("GT"); // always include GT first

    map<string, vector<string> > full_null_val; // values to set missing genotypes ('.') with
    full_null_val.insert(pair<string,vector<string> >("GT", null_GTval_vector));

    for(const auto& formats : variantFile.formatTypes){
        if( formats.first != "GT" ){
            out_formats.push_back(formats.first);
            full_null_val.insert(pair<string,vector<string> >(formats.first, null_val_vector));
        }
    }

    int last_num_alleles = -1;
    while (variantFile.getNextVariant(var)) {

        // make the expanded null/empty GL value if we need it (varies by num alleles)
        if( expand_GL && last_num_alleles != var.alleles.size() ){
            last_num_alleles = var.alleles.size();
            null_GLval_vector.back() = "0";
            for( int i=1; i<nChoosek(var.alleles.size()+ploidy-1, ploidy); ++i ){
                null_GLval_vector.back() += ",0";
            }
            full_null_val["GL"] = null_GLval_vector;
        }

        for( vector<string>::iterator sample_name=var.sampleNames.begin(); sample_name != var.sampleNames.end(); ++sample_name ){
            map<string, map<string, vector<string> > >::iterator sample_it = var.samples.find(*sample_name);
            if( sample_it == var.samples.end() ){
                // add missing ('.') samples back with null values
                var.samples.insert(pair<string, map<string, vector<string> > >(*sample_name, full_null_val));
            }else{
                for( vector<string>::iterator fit = out_formats.begin(); fit != out_formats.end(); ++fit ){
                    if( sample_it->second.find(*fit) == sample_it->second.end() ){ // missing key, should't happen but check anyway
                        vector<string>& val = null_val_vector; // fill any missing keys with null values
                    }
                    if( expand_GL ){
                        map<string, vector<string> >::iterator val_it = sample_it->second.find("GL");
                        if( val_it == sample_it->second.end() || val_it->second.back() == "." ){
                            sample_it->second["GL"] = null_GLval_vector;
                        }
                    }
                }
            }

        }

        var.format = out_formats; // set FORMAT so it is the same for all variants
        cout << var << endl;
    }
    return 0;

}
