#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
using namespace vcf;

bool samplesDiffer(vector<string>& samples, Variant& var) {

    string genotype;

    for (vector<string>::iterator s = samples.begin(); s != samples.end(); ++s) {
        string& sampleName = *s;
        map<string, map<string, vector<string> > >::iterator f = var.samples.find(sampleName);
        if (f != var.samples.end()) {
            map<string, vector<string> >& sample = f->second;
            map<string, vector<string> >::iterator gt = sample.find("GT");
            if (gt != sample.end()) {
                string& thisGenotype = gt->second.front();
                if (genotype.empty()) {
                    genotype = thisGenotype;
                } else {
                    if (genotype != thisGenotype) {
                        return true;
                    }
                }
            }
        }
    }

    return false;

}


int main(int argc, char** argv) {

    if (argc < 5) {
        cerr << "usage: " << argv[0] << " <tag> <sample> <sample> [ <sample> ... ] <vcf file>" << endl
             << "tags each record where the listed sample genotypes differ with <tag>" << endl;
        return 1;
    }

    string tag = argv[1];

    vector<string> samples;
    for (int i = 2; i < argc - 1; ++i) {
        samples.push_back(argv[i]);
    }

    string filename = argv[argc-1];

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

    // TODO check if AC is present
    // ensure that AC is listed as an info field
    string line = "##INFO=<ID=" + tag + ",Number=1,Type=Flag,Description=\"Samples";
    for (vector<string>::iterator s = samples.begin(); s != samples.end(); ++s) {
        line += " " + *s;
    }
    line += " have different genotypes\">";
    variantFile.addHeaderLine(line);

    // write the new header
    cout << variantFile.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        if (samplesDiffer(samples, var)) {
            var.infoFlags[tag] = true;
        }
        cout << var << endl;
    }

    return 0;

}

