#include "Variant.h"
#include <getopt.h>
#include "Fasta.h"
#include "gpatInfo.hpp"
#include <algorithm>
#include <list>
#include <set>

using namespace std;
using namespace vcflib;


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [<vcf file>]" << endl
         << endl
         << "options:" << endl 
         << "    -h, --help              Print this message" << endl
         << "    -v, --version           Print version" << endl
         << "    -r, --reference FILE    FASTA reference file, required with -i and -u" << endl
         << "    -w, --window-size N     Merge variants at most this many bp apart (default 30)" << endl
         << "    -o, --only-variants     Don't output the entire haplotype, just concatenate" << endl
         << "                            REF/ALT strings (delimited by \":\")" << endl
         << endl
         << "Convert genotype-based phased alleles within --window-size into haplotype alleles." << endl
         << "Will break haplotype construction when encountering non-phased genotypes on input." << endl
         << endl;
    exit(0);
}

bool isPhased(Variant& var) {
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator g = sample.find("GT");
        if (g != sample.end()) {
            string gt = g->second.front();
            if (gt.size() > 1 && gt.find("|") == string::npos) {
                return false;
            }
        }
    }
    return true;
}

int main(int argc, char** argv) {

    string vcfFileName;
    string fastaFileName;
    int windowsize = 30;
    bool onlyVariants = false;

    if (argc == 1)
        printSummary(argv);

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"window-size", required_argument, 0, 'w'},
                {"reference", required_argument, 0, 'r'},
                {"version", no_argument, 0, 'v'},
                {"only-variants", no_argument, 0, 'o'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvow:r:",
                         long_options, &option_index);
	
        if (c == -1){
            break;
        }
        switch (c) {
        case 'v':
        {
            printBasicVersion();
            exit(0);
        }  
        case 'o':
        {
            onlyVariants = true;
            break;
        }
        case 'w':
        {
            windowsize = atoi(optarg);
            break;
        }
        case 'r':
        {
            fastaFileName = string(optarg);
            break;
        }
        case 'h':
        {
            printSummary(argv);
            break;
        }
        case '?':
        {
            printSummary(argv);
            exit(1);
            break;
        }
        default:
            abort ();
        }
    }

    VariantCallFile variantFile;
    string inputFilename;
    if (optind == argc - 1) {
        inputFilename = argv[optind];
        variantFile.open(inputFilename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        cerr << "could not open VCF file" << endl;
        exit(1);
    }

    FastaReference reference;
    if (fastaFileName.empty()) {
        cerr << "a reference is required for haplotype allele generation" << endl;
        exit(1);
    }
    reference.open(fastaFileName);

    // pattern
    // when variants are within windowSize from each other, build up local haplotypes
    // establish all the haplotypes which exist within the window using genotypes+allele#+position map
    // generate a haplotype allele string for each unique haplotype
    // for completeness retain phasing information in the genotypes
    // write a new VCF record in which there are haplotype alleles and correctly described genotypes for each sample
    // if the variants are outside of the windowSize, just write out the record

    Variant var(variantFile);
    Variant outputVar(variantFile);

    cout << variantFile.header << endl;

    // get the first distances
    vector<Variant> cluster;

    while (variantFile.getNextVariant(var) || !cluster.empty()) {

        bool haplotypeCluster = false;

        if (variantFile.done()) {
            if (cluster.size() >= 1) {
                haplotypeCluster = true;
            } else {
                cout << cluster.front() << endl;
                cluster.clear();
            }
        } else if (isPhased(var)) {
            if (cluster.empty()
                || cluster.back().sequenceName == var.sequenceName
                && var.position - cluster.back().position + cluster.back().ref.size() - 1 <= windowsize) {
                cluster.push_back(var);
            } else {
                if (cluster.size() == 1) {
                    cout << cluster.front() << endl;
                    cluster.clear();
                    if (!variantFile.done()) {
                        cluster.push_back(var);
                    }
                } else {
                    haplotypeCluster = true;
                }
            }
        } else { // not phased
            if (cluster.empty()) {
                cout << var << endl;
            } else if (cluster.size() == 1) {
                cout << cluster.front() << endl;
                cluster.clear();
                cout << var << endl;
            } else {
                haplotypeCluster = true;
            }
        }

        // we need to deal with the current cluster, as our next var is outside of bounds
        // process the last cluster if it's more than 1 var
        if (haplotypeCluster) {
            /*
            cerr << "cluster: ";
            for (vector<Variant>::iterator v = cluster.begin(); v != cluster.end(); ++v) {
                cerr << " " << v->position;
            }
            cerr << endl;
            */

            // generate haplotype alleles and genotypes!
            // get the reference sequence across the haplotype in question
            string referenceHaplotype = reference.getSubSequence(cluster.front().sequenceName,
                                                                 cluster.front().position - 1,
                                                                 cluster.back().position
                                                                 + cluster.back().ref.size() - cluster.front().position);

            // establish what haplotypes there are by parsing the (phased) genotypes across the samples over these records
            map<string, vector<vector<int> > > sampleHaplotypes;
            for (vector<string>::iterator s = var.sampleNames.begin(); s != var.sampleNames.end(); ++s) {
                // build the haplotype using the genotype fields in the variant cluster
                // only build haplotypes for samples with complete information
                string& sampleName = *s;
                vector<vector<int> >& haplotypes = sampleHaplotypes[sampleName];
		
                bool completeCoverage = true;
                // ensure complete genotype coverage over the haplotype cluster
                for (vector<Variant>::iterator v = cluster.begin(); v != cluster.end(); ++v) {
                    if (v->samples.find(sampleName) == v->samples.end()
                        || v->samples[sampleName].find("GT") == v->samples[sampleName].end()) {
                        completeCoverage = false;
                        break;
                    }
                }
                if (!completeCoverage) {
                    continue; // skip samples without complete coverage
                }
		
                // what's the ploidy?
                {
                    string& gt = cluster.front().samples[sampleName]["GT"].front();
                    vector<string> gtspec = split(gt, "|");
                    for (vector<string>::iterator g = gtspec.begin(); g != gtspec.end(); ++g) {
                        vector<int> haplotype;
                        haplotypes.push_back(haplotype);
                    }
                }
		
                for (vector<Variant>::iterator v = cluster.begin(); v != cluster.end(); ++v) {
                    string& gt = v->samples[sampleName]["GT"].front();
                    vector<string> gtspec = split(gt, "|");
                    vector<string>::iterator g = gtspec.begin();
                    for (vector<vector<int> >::iterator h = haplotypes.begin(); h != haplotypes.end(); ++h, ++g) {
                        int j;
                        convert(*g, j);
                        h->push_back(j);
                    }
                }
            }

            map<const vector<int>*, vector<string> > hapToSamples;
            set<vector<int> > uniqueHaplotypes;
            for (map<string, vector<vector<int> > >::iterator hs = sampleHaplotypes.begin();
                 hs != sampleHaplotypes.end(); ++hs) {
                vector<vector<int> >& haps = hs->second;
                for (vector<vector<int> >::iterator h = haps.begin(); h != haps.end(); ++h) {
                    uniqueHaplotypes.insert(*h);
                    hapToSamples[&*uniqueHaplotypes.find(*h)].push_back(hs->first);
                }
            }
	    
            // write new haplotypes
            map<vector<int>, string> haplotypeSeqs;
            map<vector<int>, int> haplotypeIndexes;
            map<int, string> alleles;
	    
            int impossibleHaplotypes = 0;

            // always include the reference haplotype as 0
            // when we come to it in the haplotypes, we'll ignore it
            int alleleIndex = 1;
            for (set<vector<int> >::iterator u = uniqueHaplotypes.begin(); u != uniqueHaplotypes.end(); ++u) {

                /*
                for (vector<int>::const_iterator z = u->begin(); z != u->end(); ++z) {
                    cerr << *z;
                }
                cerr << endl;
                */

                string haplotype;
                if (!onlyVariants) {
                    haplotype = referenceHaplotype;
                }
                bool isreference = true;
                bool impossibleHaplotype = false;
                int referenceInsertOffset = 0;
                int j = 0; // index into variant cluster
                int lastpos = 0;
                int lastrefend = 0;
                for (vector<int>::const_iterator z = u->begin(); z != u->end(); ++z, ++j) {
                    int i = *z;
                    Variant& vartoInsert = cluster.at(j);
                    if (i == 0) {
                        if (onlyVariants) {
                            if (!haplotype.empty()) haplotype.append(":");
                            haplotype.append(vartoInsert.ref);
                        }
                    }
                    if (i != 0) {
                        isreference = false;
                        string& alternate = vartoInsert.alleles.at(i);
                        if (vartoInsert.position < lastrefend) {
                            cerr << "impossible haplotype, overlapping alleles at " << vartoInsert.sequenceName << ":" << vartoInsert.position << endl;
                            cerr << "+target " << vartoInsert.sequenceName << " " << vartoInsert.position-1 << " " << vartoInsert.position-1 + vartoInsert.ref.size() << endl;
                            cerr << "+variant " << vartoInsert.sequenceName << ":" << vartoInsert.position << ":" << alternate << endl;
                            impossibleHaplotype = true;
                            // find the impossible haplotype samples
                            cerr << "+samples ";
                            for (auto& sample : hapToSamples[&*u]) {
                                cerr << sample << " ";
                            } cerr << endl;
                            break;
                        } else {
                            //cerr << vartoInsert.position << " " << cluster.front().position + referenceInsertOffset << endl;
                            //cerr << "replacing " << vartoInsert.ref << " at " << vartoInsert.position - cluster.front().position + referenceInsertOffset << " with " << alternate << endl;
                            if (onlyVariants) {
                                if (!haplotype.empty()) haplotype.append(":");
                                haplotype.append(alternate);
                            } else {
                                haplotype.replace(vartoInsert.position - cluster.front().position + referenceInsertOffset,
                                                  vartoInsert.ref.size(), alternate);
                                if (alternate.size() != vartoInsert.ref.size()) {
                                    referenceInsertOffset += alternate.size() - vartoInsert.ref.size();
                                }
                                lastpos = vartoInsert.position;
                                lastrefend = vartoInsert.position + vartoInsert.ref.size();
                            }
                        }
                    }
                }
		
                if (impossibleHaplotype) {
                    ++impossibleHaplotypes;
                    haplotypeIndexes[*u] = -1; // indicates impossible haplotype
                    impossibleHaplotype = false;
                } else if (isreference) {
                    alleles[0] = haplotype;
                    haplotypeIndexes[*u] = 0;
                } else {
                    alleles[alleleIndex] = haplotype;
                    haplotypeIndexes[*u] = alleleIndex;
                    ++alleleIndex;
                }
                haplotypeSeqs[*u] = haplotype;
                // if there's not a reference allele, add it
                if (alleles.find(0) == alleles.end()) {
                    alleles[0] = referenceHaplotype;
                    // nb, there is no reference haplotype among
                    // the samples, so we don't have to add it to
                    // the haplotypeIndexes
                }
            }

            if (onlyVariants) {
                string newRef;
                for (vector<Variant>::iterator v = cluster.begin(); v != cluster.end(); ++v) {
                    if (!newRef.empty()) newRef.append(":");
                    newRef.append(v->ref);
                }
                outputVar.ref = newRef;
            } else {
                outputVar.ref = alleles[0];
            }
            outputVar.alt.clear();
            for (int i = 1; i < alleleIndex; ++i) {
                outputVar.alt.push_back(alleles[i]);
            }
	    
            outputVar.sequenceName = cluster.front().sequenceName;
            outputVar.position = cluster.front().position;
            outputVar.filter = ".";
            outputVar.id = ".";
            outputVar.info = cluster.front().info;
            outputVar.samples.clear();
            outputVar.format = cluster.front().format;
	    
            // now the genotypes
            for (vector<string>::iterator s = var.sampleNames.begin(); s != var.sampleNames.end(); ++s) {
                string& sampleName = *s;
                vector<string> gt;
                vector<vector<int> > & hs = sampleHaplotypes[sampleName];
                for (vector<vector<int> >::iterator h = hs.begin(); h != hs.end(); ++h) {
                    int hi = haplotypeIndexes[*h];
                    if (hi != -1) {
                        gt.push_back(convert(hi));
                    } else {
                        // nonexistent or impossible haplotype
                        gt.push_back(".");
                    }
                }
                if (gt.size() != 0) {
                    outputVar.samples[sampleName]["GT"].push_back(join(gt, "|"));
                }
            }
            if (cluster.size() - impossibleHaplotypes < 2) {
                for (vector<Variant>::iterator v = cluster.begin(); v != cluster.end(); ++v) {
                    cout << *v << endl;
                }
            } else {
                if (!outputVar.alt.empty()) {
                    cout << outputVar << endl;
                } else {
                    cerr << "no alternate alleles remain at " << outputVar.sequenceName << ":" << outputVar.position << " after haplotype validation" << endl;
                }
            }
            cluster.clear();
            if (!variantFile.done()) cluster.push_back(var);
        }
    }

    exit(0);  // why?
    return 0;

}

