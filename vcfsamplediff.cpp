#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>
#include <getopt.h>

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
             << "tags each record where the listed sample genotypes differ with <tag>" << endl
	     << "The first sample is assumed to be germline, the second somatic." << endl
	     << "Each record is tagged with <tag>={germline,somatic,loh} to specify the type of" << endl
	     << "variant given the genotype difference between the two samples." << endl;
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

    assert(samples.size() == 2);

    Variant var(variantFile);

    // TODO check if AC is present
    // ensure that AC is listed as an info field
    string line = "##INFO=<ID=" + tag + ",Number=1,Type=String,Description=\"Samples";
    for (vector<string>::iterator s = samples.begin(); s != samples.end(); ++s) {
        line += " " + *s;
    }
    line += " have different genotypes\">";
    variantFile.addHeaderLine(line);

    // write the new header
    cout << variantFile.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
	if (var.samples.find(samples.front()) != var.samples.end()
	    && var.samples.find(samples.back()) != var.samples.end()) {
	    map<int, int> gtGermline = decomposeGenotype(var.samples[samples.front()]["GT"].front());
	    map<int, int> gtSomatic  = decomposeGenotype(var.samples[samples.back()]["GT"].front());
	    if (gtGermline == gtSomatic) {
		var.info[tag].push_back("germline");
	    } else {
		if (isHet(gtGermline) && isHom(gtSomatic)) {
		    var.info[tag].push_back("loh");
		} else if (isHom(gtGermline) && isHet(gtSomatic)) {
		    var.info[tag].push_back("somatic");
		}
	    }
	}
        cout << var << endl;
    }

    return 0;

}

