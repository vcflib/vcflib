#include "Variant.h"
#include "split.h"
#include <string>
#include <set>
#include <iostream>

using namespace std;
using namespace vcf;


int main(int argc, char** argv) {

    if (argc != 3) {
        cerr << "usage: " << argv[0] << " <vcf file>" << endl
             << "report the allele freqs and call rates for each variant in the vcf file" << endl;
        return 1;
    }

    string vcffile    = argv[1];
    string samplefile  = argv[2];

    VariantCallFile variantFile;

    if (vcffile == "-") {
        variantFile.open(std::cin);
    } else {
        variantFile.open(vcffile);
    }

    if (!variantFile.is_open()) {
        return 1;
    }
    
    
    ifstream samples(samplefile.c_str(), ios::in);
    string sample;
    string ethnicity;
    map<string, string> sample_to_cohort;
    map<string, string> sample_to_ethnicity;
    set<string> cohorts;
    set<string> ethnicities;
    vector<string> sample_pieces;
    while (samples >> sample >> ethnicity) {
        sample_pieces.clear();
        split(sample, ":", sample_pieces);
        // cohort is the first piece of the split string, e.g., HeartGO:2333
        string cohort = sample_pieces[0];
        sample_to_cohort[sample]    = cohort;
        sample_to_ethnicity[sample] = ethnicity;
        cohorts.insert(cohort);
        ethnicities.insert(ethnicity);
    }
    
    // MAKE A HEADER
    cout << "chrom\tstart\tend\tref\talt\tid\tfilter\tquality\tinfo\tgene\timpact\tnum_samples\tnum_valid_genotypes\ttotal_called_alleles\ttotal_alt_alleles\t";
    set<string>::iterator c     = cohorts.begin();
    for (; c != cohorts.end(); ++c) {
        cout << *c << "-alt-called\t" << *c << "-tot-called\t";
    }
    cout << "\t";
    set<string>::iterator e     = ethnicities.begin();
    for (; e != ethnicities.end(); ++e) {
        cout << *e << "-alt-called\t" << *e << "-tot-called\t";
    }
    cout << endl;
    
    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        map<string, map<string, vector<string> > >::iterator s     = var.samples.begin(); 
        map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();
        
        cout << var.sequenceName << "\t"
             << var.position-1     << "\t"
             << var.position     << "\t"
             << var.ref          << "\t";
        var.printAlt(cout);     cout << "\t"; 
        cout << var.filter << "\t";
        cout << var.id << "\t";
        cout << var.quality << "\t";

        // store the set of genes, transcripts, and impacts;
        
        map<string, vector<string> >::iterator i     = var.info.begin();
        map<string, vector<string> >::iterator iEnd  = var.info.end();
        for (; i != iEnd; ++i) {
            if (i->first == "DP" ||
                i->first == "NS" ||
                i->first == "AM" ||
                i->first == "AL" ||
                i->first == "AF" ||
                i->first == "AB" ||
                i->first == "STR" ||
                i->first == "STZ" ||
                i->first == "FIC" ||
                i->first == "ANNO") 
            {
                cout << i->first << ":";
                vector<string> values = i->second;
                for (size_t v = 0; v < values.size(); ++v) {
                    cout << values[v];
                    cout << ";";
                }
            }
        }
        // extract the gene and impact from the ANNO info tag
        string gene = "none";
        string impact = "none";
        i     = var.info.begin();
        for (; i != iEnd; ++i) {
            if (i->first == "ANNO") 
            {
                vector<string> gene_anno = split(i->second[0], ":");
                gene        = gene_anno[1];
                impact      = gene_anno[0];
            }
        }
        
        cout << "\t" << gene << "\t" << impact;
        
        // count the number of successful (i.e., not ./.) genotypes
        int valid_genotypes = 0;
        int alt_alleles = 0;
        int called_alleles = 0;
        
        map<string, int> cohort_tot_called_alleles;
        map<string, int> cohort_tot_alt_alleles;
        map<string, int> eth_tot_called_alleles;
        map<string, int> eth_tot_alt_alleles;
        for (; s != sEnd; ++s) {
            
            string sample = s->first;
            string cohort    = sample_to_cohort[sample];
            string ethnicity = sample_to_ethnicity[sample];
            
            map<string, vector<string> >& sample_info = s->second;
            string genotype       = sample_info["GT"].front();
            vector<string> alleles = split(genotype, "|/");
            //cout << genotype << "=";
            for (size_t i = 0; i < alleles.size(); ++i) {
                //cout << alleles[i] << ":";
                if (alleles[i] != "0" && alleles[i] != ".") {
                    alt_alleles++;
                    cohort_tot_alt_alleles[cohort]++;
                    eth_tot_alt_alleles[ethnicity]++;
                    //cout << "alt";
                }
                if (alleles[i] != ".") {
                    called_alleles++;
                    cohort_tot_called_alleles[cohort]++;
                    eth_tot_called_alleles[ethnicity]++;
                    //cout << "tot";
                }
                //cout << ";";
            }
            if (genotype != "./.")
                valid_genotypes++;
        }
        cout << "\t" 
             << var.samples.size() << "\t"
             << valid_genotypes << "\t"
             << called_alleles << "\t" 
             << alt_alleles << "\t";
        c     = cohorts.begin();
        for (; c != cohorts.end(); ++c) {
            cout << cohort_tot_alt_alleles[*c] << "\t" << cohort_tot_called_alleles[*c] << "\t";
        }
        e     = ethnicities.begin();
        for (; e != ethnicities.end(); ++e) {
            cout << eth_tot_alt_alleles[*e] << "\t" << eth_tot_called_alleles[*e] << "\t";
        }
        cout << endl;
    }
    return 0;

}

