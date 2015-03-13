#include "Variant.h"
#include "convert.h"

using namespace std;
using namespace vcflib;


double convertStrDbl(const string& s) {
    double r;
    convert(s, r);
    return r;
}

int main(int argc, char** argv) {

    int maxAlleles = 2;

    VariantCallFile variantFile;

    if (argc > 1) {
        string filename = argv[1];
        if (filename == "--help" || filename == "-h") {
            cerr << "usage: vcfflatten [file]" << endl
                 << endl
                 << "Removes multi-allelic sites by picking the most common alternate.  Requires" << endl
                 << "allele frequency specification 'AF' and use of 'G' and 'A' to specify the" << endl
                 << "fields which vary according to the Allele or Genotype. VCF file may be" << endl
                 << "specified on the command line or piped as stdin." << endl;
            exit(1);
        }
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
        // count the number of alternates
        // if we have more than N, strip the lowest-frequency ones
        if (var.alleles.size() > maxAlleles) {

            multimap<double, string> alleleFrequencies;

            vector<string>& freqsstr = var.info["AF"];
            vector<double> freqs;
            freqs.resize(freqsstr.size());
            transform(freqsstr.begin(), freqsstr.end(), freqs.begin(), convertStrDbl);

            vector<double>::iterator f = freqs.begin();
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a, ++f) {
                alleleFrequencies.insert(pair<double, string>(*f, *a));
            }

            // pick the highest frequency alternate
            string bestalt = alleleFrequencies.rbegin()->second;
            // and get its index
            int bestaltIndex = var.getAltAlleleIndex(bestalt);
            int bestaltGenotypeIndex = bestaltIndex + 1; // per VCF spec

            // keep the RR, RA, and AA alleles for this alternate
            // generate the genotype index table for this variant
            map<pair<int, int>, int> genotypeIndexes = var.getGenotypeIndexesDiploid();

            // now get the genotype indexes we want to keep
            vector<int> alleleIndexes;
            alleleIndexes.push_back(0); 
            alleleIndexes.push_back(bestaltGenotypeIndex);

            // add the reference allele index for generating genotype indexes
            int ploidy = 2;
            vector<vector<int> > genotypesToKeep = multichoose(ploidy, alleleIndexes);
            map<int, bool> genotypeIndexesToKeep;
            for (vector<vector<int> >::iterator k = genotypesToKeep.begin(); k != genotypesToKeep.end(); ++k) {
                pair<int, int> genotype = make_pair(k->front(), k->back()); // vectors are guaranteed to be diploid per multichoose
                genotypeIndexesToKeep[genotypeIndexes[genotype]] = true;
            }
            // we are diploid, so there should be exactly 3 genotypes
            assert(genotypeIndexesToKeep.size() == 3);

            // get the fields which have genotype order "G"
            // for all the infocounts
            // find the ones which are == GENOTYPE_NUMBER or ALLELE_NUMBER
            //     and fix em up
            for (map<string, int>::iterator c = variantFile.infoCounts.begin(); c != variantFile.infoCounts.end(); ++c) {
                int count = c->second;
                if (count == GENOTYPE_NUMBER) {
                    string key = c->first;
                    map<string, vector<string> >::iterator v = var.info.find(key);
                    if (v != var.info.end()) {
                        vector<string>& vals = v->second;
                        vector<string> tokeep;
                        int i = 0;
                        for (vector<string>::iterator g = vals.begin(); g != vals.end(); ++g, ++i) {
                            if (genotypeIndexesToKeep.find(i) != genotypeIndexesToKeep.end()) {
                                tokeep.push_back(*g);
                            }
                        }
                        vals = tokeep;
                    }
                } else if (count == ALLELE_NUMBER) {
                    string key = c->first;
                    map<string, vector<string> >::iterator v = var.info.find(key);
                    if (v != var.info.end()) {
                        vector<string>& vals = v->second;
                        vector<string> tokeep;
                        int i = 0;
                        for (vector<string>::iterator a = vals.begin(); a != vals.end(); ++a, ++i) {
                            if (i == bestaltIndex) {
                                tokeep.push_back(*a);
                            }
                        }
                        vals = tokeep;
                    }
                }
            }
            //
            // for all the formatcounts
            // find the ones which are == GENOTYPE_NUMBER or ALLELE_NUMBER
            //     for each sample, remove the new irrelevant values

            // for each sample
            //   remove info fields which now refer to nothing
            for (map<string, int>::iterator c = variantFile.formatCounts.begin(); c != variantFile.formatCounts.end(); ++c) {
                int count = c->second;
                if (count == GENOTYPE_NUMBER) {
                    string key = c->first;
                    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
                        map<string, vector<string> >& sample = s->second;
                        map<string, vector<string> >::iterator v = sample.find(key);
                        if (v != sample.end()) {
                            vector<string>& vals = v->second;
                            vector<string> tokeep;
                            int i = 0;
                            for (vector<string>::iterator g = vals.begin(); g != vals.end(); ++g, ++i) {
                                if (genotypeIndexesToKeep.find(i) != genotypeIndexesToKeep.end()) {
                                    tokeep.push_back(*g);
                                }
                            }
                            vals = tokeep;
                        }
                    }
                } else if (count == ALLELE_NUMBER) {
                    string key = c->first;
                    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
                        map<string, vector<string> >& sample = s->second;
                        map<string, vector<string> >::iterator v = sample.find(key);
                        if (v != sample.end()) {
                            vector<string>& vals = v->second;
                            vector<string> tokeep;
                            int i = 0;
                            for (vector<string>::iterator a = vals.begin(); a != vals.end(); ++a, ++i) {
                                if (i == bestaltIndex) {
                                    tokeep.push_back(*a);
                                }
                            }
                            vals = tokeep;
                        }
                    }
                }
            }

            var.alt.clear();
            var.alt.push_back(bestalt);

        }
        cout << var << endl;
    }

    return 0;

}

