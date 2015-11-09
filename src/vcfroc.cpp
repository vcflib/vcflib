#include "Variant.h"
#include "BedReader.h"
#include "IntervalTree.h"
#include <getopt.h>
#include "Fasta.h"
#include <algorithm>
#include <list>
#include <set>

using namespace std;
using namespace vcflib;


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [<vcf file>]" << endl
         << endl
         << "options:" << endl 
         << "    -t, --truth-vcf FILE      use this VCF as ground truth for ROC generation" << endl
         << "    -w, --window-size N       compare records up to this many bp away (default 30)" << endl
         << "    -c, --complex             directly compare complex alleles, don't parse into primitives" << endl
         << "    -r, --reference FILE      FASTA reference file" << endl
         << endl
         << "Generates a pseudo-ROC curve using sensitivity and specificity estimated against" << endl
         << "a putative truth set.  Thresholding is provided by successive QUAL cutoffs." << endl;
    exit(0);
}

void buildVariantIntervalTree(VariantCallFile& variantFile,
                              map<string, IntervalTree<Variant*> >& variantIntervals,
                              list<Variant>& variants) {

    map<string, vector<Interval<Variant*> > > rawVariantIntervals;
    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
        long int left = var.position;
        long int right = left + var.ref.size(); // this should be 1-past the end
        variants.push_back(var);
        Variant* v = &variants.back();
        rawVariantIntervals[var.sequenceName].push_back(Interval<Variant*>(left, right, v));
    }
	
    for (map<string, vector<Interval<Variant*> > >::iterator j = rawVariantIntervals.begin(); j != rawVariantIntervals.end(); ++j) {
        variantIntervals[j->first] = IntervalTree<Variant*>(j->second);
    }
}


void intersectVariant(Variant& var,
                      map<string, IntervalTree<Variant*> >& variantIntervals,
                      vector<string*>& commonAlleles,
                      vector<string*>& uniqueAlleles,
                      FastaReference& reference,
                      int windowsize = 50) {

    vector<Interval<Variant*> > results;

    variantIntervals[var.sequenceName].findContained(var.position - windowsize, var.position + var.ref.size() + windowsize, results);

    vector<Variant*> overlapping;

    for (vector<Interval<Variant*> >::iterator r = results.begin(); r != results.end(); ++r) {
        overlapping.push_back(r->value);
    }


    if (overlapping.empty()) {
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            uniqueAlleles.push_back(&*a);
        }
    } else {

        // get the min and max of the overlaps

        int haplotypeStart = var.position;
        int haplotypeEnd = var.position + var.ref.size();

        for (vector<Variant*>::iterator v = overlapping.begin(); v != overlapping.end(); ++v) {
            haplotypeStart = min((*v)->position, (long int) haplotypeStart);
            haplotypeEnd = max((*v)->position + (*v)->ref.size(), (long unsigned int) haplotypeEnd);
        }

        // for everything overlapping and the current variant, construct the local haplotype within the bounds
        // if there is an exact match, the allele in the current VCF does intersect

        string referenceHaplotype = reference.getSubSequence(var.sequenceName, haplotypeStart - 1, haplotypeEnd - haplotypeStart);
        map<string, vector<pair<Variant*, int> > > haplotypes; // map to variant and alt index

        for (vector<Variant*>::iterator v = overlapping.begin(); v != overlapping.end(); ++v) {
            Variant& variant = **v;
            int altindex = 0;
            for (vector<string>::iterator a = variant.alt.begin(); a != variant.alt.end(); ++a, ++altindex) {
                string haplotype = referenceHaplotype;
                // get the relative start and end coordinates for the variant alternate allele
                int relativeStart = variant.position - haplotypeStart;
                haplotype.replace(relativeStart, variant.ref.size(), *a);
                haplotypes[haplotype].push_back(make_pair(*v, altindex));
            }
        }


        // determine the non-intersecting alts
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            string haplotype = referenceHaplotype;
            int relativeStart = var.position - haplotypeStart;
            haplotype.replace(relativeStart, var.ref.size(), *a);
            map<string, vector<pair<Variant*, int> > >::iterator h = haplotypes.find(haplotype);
            if (h == haplotypes.end()) {
                uniqueAlleles.push_back(&*a);
            } else {
                commonAlleles.push_back(&*a);
            }
        }

    }
}


int main(int argc, char** argv) {

    string truthVcfFileName;
    string fastaFileName;
    bool complex = false;
    int windowsize = 30;

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
                {"complex", required_argument, 0, 'c'},
                {"truth-vcf", required_argument, 0, 't'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hcw:r:t:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

	    case 'w':
            windowsize = atoi(optarg);
            break;

	    case 'r':
            fastaFileName = string(optarg);
            break;

	    case 't':
	        truthVcfFileName = optarg;
            break;

        case 'c':
            complex = true;
            break;

        case 'h':
            printSummary(argv);
            break;

        case '?':
            printSummary(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    VariantCallFile variantFile;
    bool usingstdin = false;
    string inputFilename;
    if (optind == argc - 1) {
        inputFilename = argv[optind];
        variantFile.open(inputFilename);
    } else {
        variantFile.open(std::cin);
        usingstdin = true;
    }

    if (!variantFile.is_open()) {
        cerr << "could not open VCF file" << endl;
        exit(1);
    }

    VariantCallFile truthVariantFile;
    if (!truthVcfFileName.empty()) {
        if (truthVcfFileName == "-") {
            if (usingstdin) {
                cerr << "cannot open both VCF file streams from stdin" << endl;
                exit(1);
            } else {
                truthVariantFile.open(std::cin);
            }
        } else {
            truthVariantFile.open(truthVcfFileName);
        }
        if (!truthVariantFile.is_open()) {
            cerr << "could not open VCF file " << truthVcfFileName << endl;
            exit(1);
        }
    }

    FastaReference reference;
    if (fastaFileName.empty()) {
        cerr << "a reference is required for the haplotype-based intersection used by vcfroc" << endl;
        exit(1);
    }
    reference.open(fastaFileName);

    // read the VCF file for union or intersection into an interval tree
    // indexed using some proximity window

    map<string, IntervalTree<Variant*> > truthVariantIntervals;
    list<Variant> truthVariants;
    buildVariantIntervalTree(truthVariantFile, truthVariantIntervals, truthVariants);

    map<string, IntervalTree<Variant*> > testVariantIntervals;
    list<Variant> testVariants;
    buildVariantIntervalTree(variantFile, testVariantIntervals, testVariants);

    map<long double, vector<VariantAllele*> > falseNegativeAllelesAtCutoff;  // false negative after this cutoff
    map<long double, vector<VariantAllele*> > falsePositiveAllelesAtCutoff;  // false positive until this cutoff
    list<VariantAllele*> allFalsePositiveAlleles;
    map<long double, vector<VariantAllele*> > allelesAtCutoff;
    //map<long double, vector<VariantAllele*> > totalAllelesAtCutoff;
    map<Variant*, map<string, vector<VariantAllele> > > parsedAlleles;
    map<long double, vector<Variant*> > callsByCutoff;

    // replicate this method, where Q is for each unique Q in the set
    //vcfintersect -r $reference -v -i $results.$Q.vcf $answers_primitives | vcfstats >false_negatives.$Q.stats
    //vcfintersect -r $reference -v -i $answers_primitives $results.$Q.vcf | vcfstats >false_positives.$Q.stats

    for (list<Variant>::iterator v = testVariants.begin(); v != testVariants.end(); ++v) {
        // TODO allow different cutoff sources
        callsByCutoff[v->quality].push_back(&*v);
    }

    // add false negatives at any cutoff
    for (list<Variant>::iterator v = truthVariants.begin(); v != truthVariants.end(); ++v) {
        Variant& variant = *v;
        vector<string*> commonAlleles;
        vector<string*> uniqueAlleles;
        intersectVariant(variant, testVariantIntervals,
                         commonAlleles, uniqueAlleles, reference);
        if (complex) {
            parsedAlleles[&*v] = variant.flatAlternates();
        } else {
            parsedAlleles[&*v] = variant.parsedAlternates();
        }
        // unique alleles are false negatives regardless of cutoff
        for (vector<string*>::iterator a = uniqueAlleles.begin(); a != uniqueAlleles.end(); ++a) {
            vector<VariantAllele>& alleles = parsedAlleles[&*v][**a];
            for (vector<VariantAllele>::iterator va = alleles.begin(); va != alleles.end(); ++va) {
                if (va->ref != va->alt) { 		// use only non-reference alleles
                    // false negatives at threshold 0 XXX --- may not apply if threshold is generalized
                    falseNegativeAllelesAtCutoff[-1].push_back(&*va);
                }
            }
        }
    }

    for (map<long double, vector<Variant*> >::iterator q = callsByCutoff.begin(); q != callsByCutoff.end(); ++q) {
        long double threshold = q->first;
        vector<Variant*>& variants = q->second;
        for (vector<Variant*>::iterator v = variants.begin(); v != variants.end(); ++v) {
            Variant& variant = **v;
            vector<string*> commonAlleles;
            vector<string*> uniqueAlleles;
            intersectVariant(variant, truthVariantIntervals,
                             commonAlleles, uniqueAlleles, reference);
            if (complex) {
                parsedAlleles[*v] = variant.flatAlternates();
            } else {
                parsedAlleles[*v] = variant.parsedAlternates();
            }

            map<string, vector<VariantAllele> >& parsedAlts = parsedAlleles[*v];
            // push VariantAllele*'s into the FN and FP alleles at cutoff vectors
            for (vector<string*>::iterator a = commonAlleles.begin(); a != commonAlleles.end(); ++a) {
                vector<VariantAllele>& alleles = parsedAlleles[*v][**a];
                for (vector<VariantAllele>::iterator va = alleles.begin(); va != alleles.end(); ++va) {
                    if (va->ref != va->alt) { 		// use only non-reference alleles
                        allelesAtCutoff[threshold].push_back(&*va);
                        falseNegativeAllelesAtCutoff[threshold].push_back(&*va);
                    }
                }
            }
            for (vector<string*>::iterator a = uniqueAlleles.begin(); a != uniqueAlleles.end(); ++a) {
                vector<VariantAllele>& alleles = parsedAlts[**a];
                for (vector<VariantAllele>::iterator va = alleles.begin(); va != alleles.end(); ++va) {
                    if (va->ref != va->alt) { 		// use only non-reference alleles
                        allelesAtCutoff[threshold].push_back(&*va);
                        allFalsePositiveAlleles.push_back(&*va);
                        falsePositiveAllelesAtCutoff[threshold].push_back(&*va);
                    }
                }
            }
        }
    }


    // output results
    int totalSNPs = 0;
    int falsePositiveSNPs = 0;
    int falseNegativeSNPs = 0;
    int totalIndels = 0;
    int falsePositiveIndels = 0;
    int falseNegativeIndels = 0;
    int totalComplex = 0;
    int falsePositiveComplex = 0;
    int falseNegativeComplex = 0;

    // write header
    
    cout << "threshold" << "\t"
         << "num_snps" << "\t"
         << "false_positive_snps" << "\t"
         << "false_negative_snps" << "\t"
         << "num_indels" << "\t"
         << "false_positive_indels" << "\t"
         << "false_negative_indels" << "\t"
         << "num_complex" << "\t"
         << "false_positive_complex" << "\t"
         << "false_negative_complex" << endl;

    // count total alleles in set
    for (map<long double, vector<VariantAllele*> >::iterator a = allelesAtCutoff.begin(); a != allelesAtCutoff.end(); ++a) {
        vector<VariantAllele*>& alleles = a->second;
        for (vector<VariantAllele*>::iterator va = alleles.begin(); va != alleles.end(); ++va) {
            VariantAllele& allele = **va;
            if (allele.ref.size() == 1 && allele.ref.size() == allele.alt.size()) {
                ++totalSNPs;
            } else if (allele.ref.size() != allele.alt.size()) {
                if (allele.ref.size() == 1 || allele.alt.size() == 1) {
                    ++totalIndels;
                } else {
                    ++totalComplex;
                }
            } else {
                ++totalComplex;
            }
        }
    }

    // tally total false positives
    for (list<VariantAllele*>::iterator va = allFalsePositiveAlleles.begin(); va != allFalsePositiveAlleles.end(); ++va) {
        VariantAllele& allele = **va;
        if (allele.ref.size() == 1 && allele.ref.size() == allele.alt.size()) {
            ++falsePositiveSNPs;
        } else if (allele.ref.size() != allele.alt.size()) {
            if (allele.ref.size() == 1 || allele.alt.size() == 1) {
                ++falsePositiveIndels;
            } else {
                ++falsePositiveComplex;
            }
        } else {
            ++falsePositiveComplex;
        }
    }

    // get categorical false negatives
    vector<VariantAllele*>& categoricalFalseNegatives = falseNegativeAllelesAtCutoff[-1];
    for (vector<VariantAllele*>::iterator va = categoricalFalseNegatives.begin(); va != categoricalFalseNegatives.end(); ++va) {
        VariantAllele& allele = **va;
        if (allele.ref.size() == 1 && allele.ref.size() == allele.alt.size()) {
            assert(allele.ref.size() == 1);
            ++falseNegativeSNPs;
        } else if (allele.ref.size() != allele.alt.size()) {
            if (allele.ref.size() == 1 || allele.alt.size() == 1) {
                ++falseNegativeIndels;
            } else {
                ++falseNegativeComplex;
            }
        } else {
            ++falseNegativeComplex;
        }
    }
    cout << -1 << "\t"
         << totalSNPs << "\t"
         << falsePositiveSNPs << "\t"
         << falseNegativeSNPs << "\t"
         << totalIndels << "\t"
         << falsePositiveIndels << "\t"
         << falseNegativeIndels << "\t"
         << totalComplex << "\t"
         << falsePositiveComplex << "\t"
         << falseNegativeComplex << endl;

    for (map<long double, vector<VariantAllele*> >::iterator a = allelesAtCutoff.begin(); a != allelesAtCutoff.end(); ++a) {
        vector<VariantAllele*>& alleles = a->second;
        long double threshold = a->first;
        for (vector<VariantAllele*>::iterator va = alleles.begin(); va != alleles.end(); ++va) {
            VariantAllele& allele = **va;
            if (allele.ref.size() == 1 && allele.ref.size() == allele.alt.size()) {
                assert(allele.ref.size() == 1);
                --totalSNPs;
            } else if (allele.ref.size() != allele.alt.size()) {
                if (allele.ref.size() == 1 || allele.alt.size() == 1) {
                    --totalIndels;
                } else {
                    --totalComplex;
                }
            } else {
                --totalComplex;
            }   
        }
        vector<VariantAllele*>& falseNegatives = falseNegativeAllelesAtCutoff[threshold];
        for (vector<VariantAllele*>::iterator va = falseNegatives.begin(); va != falseNegatives.end(); ++va) {
            VariantAllele& allele = **va;
            if (allele.ref.size() == 1 && allele.ref.size() == allele.alt.size()) {
                assert(allele.ref.size() == 1);
                ++falseNegativeSNPs;
            } else if (allele.ref.size() != allele.alt.size()) {
                if (allele.ref.size() == 1 || allele.alt.size() == 1) {
                    ++falseNegativeIndels;
                } else {
                    ++falseNegativeComplex;
                }
            } else {
                ++falseNegativeComplex;
            }
        }
        vector<VariantAllele*>& falsePositives = falsePositiveAllelesAtCutoff[threshold];
        for (vector<VariantAllele*>::iterator va = falsePositives.begin(); va != falsePositives.end(); ++va) {
            VariantAllele& allele = **va;
            if (allele.ref.size() == 1 && allele.ref.size() == allele.alt.size()) {
                assert(allele.ref.size() == 1);
                --falsePositiveSNPs;
            } else if (allele.ref.size() != allele.alt.size()) {
                if (allele.ref.size() == 1 || allele.alt.size() == 1) {
                    --falsePositiveIndels;
                } else {
                    --falsePositiveComplex;
                }
            } else {
                --falsePositiveComplex;
            }
        }
        cout << threshold << "\t"
             << totalSNPs << "\t"
             << falsePositiveSNPs << "\t"
             << falseNegativeSNPs << "\t"
             << totalIndels << "\t"
             << falsePositiveIndels << "\t"
             << falseNegativeIndels << "\t"
             << totalComplex << "\t"
             << falsePositiveComplex << "\t"
             << falseNegativeComplex << endl;

    }
    
    exit(0);  // why?
    return 0;

}

