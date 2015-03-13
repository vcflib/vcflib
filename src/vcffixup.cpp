#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
using namespace vcflib;

int countAlts(Variant& var, int alleleIndex) {
    int alts = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        if (gt != sample.end()) {
            map<int, int> genotype = decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
                if (g->first == alleleIndex) {
                    alts += g->second;
                }
            }
        }
    }
    return alts;
}

int countAlleles(Variant& var) {
    int alleles = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        if (gt != sample.end()) {
            map<int, int> genotype = decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
		if (g->first != NULL_ALLELE) {
		    alleles += g->second;
		}
            }
        }
    }
    return alleles;
}

int main(int argc, char** argv) {

  if (argc == 1 || ((argc > 1) && strcmp(argv[1], "-h") == 0) || strcmp(argv[1], "--help") == 0) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "outputs a VCF stream where AC and NS have been generated for each record using sample genotypes" << endl;
        return 1;
    }

    VariantCallFile variantFile;
    if (argc == 1 || ((argc == 2) && strcmp(argv[1], "-") == 0)) {
        variantFile.open(std::cin);
        if (!variantFile.is_open()) {
            cerr << "vcffixup: could not open stdin" << endl;
            return 1;
        }
    } else {
        string filename = argv[1];
        variantFile.open(filename);
        if (!variantFile.is_open()) {
            cerr << "vcffixup: could not open " << filename << endl;
            return 1;
        }
    }

    Variant var(variantFile);

    // remove header lines we're going to add
    variantFile.removeInfoHeaderLine("AC");
    variantFile.removeInfoHeaderLine("AF");
    variantFile.removeInfoHeaderLine("NS");
    variantFile.removeInfoHeaderLine("AN");

    // and add them back, so as not to duplicate them if they are already there
    variantFile.addHeaderLine("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">");
    variantFile.addHeaderLine("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">");
    variantFile.addHeaderLine("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">");
    variantFile.addHeaderLine("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");

    // write the new header
    cout << variantFile.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        stringstream ns;
        ns << var.samples.size();
        var.info["NS"].clear();
        var.info["NS"].push_back(ns.str());

        var.info["AC"].clear();
        var.info["AF"].clear();
        var.info["AN"].clear();

        int allelecount = countAlleles(var);
        stringstream an;
        an << allelecount;
        var.info["AN"].push_back(an.str());

        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            string& allele = *a;
            int altcount = countAlts(var, var.getAltAlleleIndex(allele) + 1);
            stringstream ac;
            ac << altcount;
            var.info["AC"].push_back(ac.str());
            stringstream af;
            af << (double) altcount / (double) allelecount;
            var.info["AF"].push_back(af.str());
        }
        cout << var << endl;
    }

    return 0;

}

