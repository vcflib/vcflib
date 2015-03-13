#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>
#include <getopt.h>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

    if (argc > 1 && (argv[1] == "-h" || argv[1] == "--help")) {
        cerr << "usage: " << argv[0] << " <[input file] >[output vcf]" << endl
             << "Adds summary statistics to each record summarizing qualities reported in" << endl
             << "called genotypes.  Uses:" << endl
             << "RO (reference observation count), QR (quality sum reference observations)" << endl
             << "AO (alternate observation count), QA (quality sum alternate observations)" << endl;
        return 1;
    }

    VariantCallFile variantFile;
    if (argc == 1) {
        variantFile.open(std::cin);
    } else {
        string filename = argv[argc-1];
        variantFile.open(filename);
        if (!variantFile.is_open()) {
            cerr << "could not open " << filename << endl;
            return 1;
        }
    }

    Variant var(variantFile);

    variantFile.removeInfoHeaderLine("AQR");
    variantFile.addHeaderLine("##INFO=<ID=AQR,Number=1,Type=Float,Description=\"Mean reference observation quality calculated by RO and QR in called samples.\">");
    variantFile.removeInfoHeaderLine("AQA");
    variantFile.addHeaderLine("##INFO=<ID=AQA,Number=A,Type=Float,Description=\"Mean alternate observation quality calculated by AO and QA in called samples.\">");
    variantFile.removeInfoHeaderLine("QR");
    variantFile.addHeaderLine("##INFO=<ID=QR,Number=1,Type=Float,Description=\"Quality sum of reference observations calculated by QR in called samples.\">");
    variantFile.removeInfoHeaderLine("QA");
    variantFile.addHeaderLine("##INFO=<ID=QA,Number=A,Type=Float,Description=\"Quality sum of alternate observations calculated by QA in called samples.\">");
    variantFile.removeInfoHeaderLine("RQA");
    variantFile.addHeaderLine("##INFO=<ID=RQA,Number=A,Type=Float,Description=\"Ratio of mean alternate observation quality to mean reference observation quality (MQA/MQR).\">");

    // write the new header
    cout << variantFile.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        int refobs = 0;
        int refqual = 0;
        vector<int> altobs(var.alt.size(), 0);
        vector<int> altqual(var.alt.size(), 0);
        for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin();
             s != var.samples.end(); ++s) {
            map<string, vector<string> >& sample = s->second;
            int x;
            if (sample.find("RO") != sample.end()) {
                convert(sample["RO"].front(), x);
                refobs += x;
                if (sample.find("QR") != sample.end()) {
                    convert(sample["QR"].front(), x);
                    refqual += x;
                }
            }
            if (sample.find("AO") != sample.end()) {
                vector<string>& aos = sample["AO"];
                for (int i = 0; i != var.alt.size(); ++i) {
                    convert(aos[i], x);
                    altobs[i] += x;
                }
                if (sample.find("QA") != sample.end()) {
                    vector<string>& qas = sample["QA"];
                    for (int i = 0; i != var.alt.size(); ++i) {
                        convert(qas[i], x);
                        altqual[i] += x;
                    }
                }
            }
        }
        var.info["QR"].push_back(convert(refqual));
        if (refobs == 0 || refqual == 0) {
            var.info["AQR"].push_back(convert(0));
        } else {
            var.info["AQR"].push_back(convert((double)refqual/(double)refobs));
        }

        for (int i = 0; i != var.alt.size(); ++i) {
            var.info["QA"].push_back(convert(altqual[i]));
            var.info["AQA"].push_back(convert((double)altqual[i]/(double)altobs[i]));
            if (refobs == 0 || refqual == 0) {
                var.info["RQA"].push_back(convert(1));
            } else {
                var.info["RQA"].push_back(convert(((double)altqual[i]/(double)altobs[i]) / 
                                                  ((double)refqual/(double)refobs)));
            }
        }
        cout << var << endl;
    }

    return 0;

}

