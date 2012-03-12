#include "Variant.h"
#include "BedReader.h"
#include "intervaltree/IntervalTree.h"
#include <getopt.h>
#include "fastahack/Fasta.h"
#include <algorithm>
#include <list>
#include <set>

using namespace std;
using namespace vcf;


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [<vcf file>]" << endl
         << endl
         << "options:" << endl 
         << "    -b, --bed FILE            use intervals provided by this BED file" << endl
         << "    -v, --invert              invert the selection, printing only records which would" << endl
         << "                                not have been printed out" << endl
	 << "    -i, --intersect-vcf FILE  use this VCF for set intersection generation" << endl
	 << "    -u, --union-vcf FILE      use this VCF for set union generation" << endl
	 << "    -w, --window-size N       compare records up to this many bp away (default 30)" << endl
	 << "    -r, --reference FILE      FASTA reference file, required with -i and -u" << endl
	 << "    -l, --loci                output whole loci when one alternate allele matches" << endl
         << endl
	 << "For bed-vcf intersection, alleles which fall into the targets are retained." << endl
	 << endl
	 << "For vcf-vcf intersection and union, unify on equivalent alleles within window-size bp" << endl
	 << "as determined by haplotype comparison alleles." << endl;
	//<< "Intersect the records in the VCF file with targets provided in a BED file." << endl
	//<< "Intersections are done on the reference sequences in the VCF file." << endl
	//<< "If no VCF filename is specified on the command line (last argument) the VCF" << endl
	//<< "read from stdin." << endl;
    exit(0);
}

int main(int argc, char** argv) {

    string bedFileName;
    string vcfFileName;
    string fastaFileName;
    bool intersecting = false;
    bool unioning = false;
    bool invert = false;
    bool contained = true;
    bool overlapping = false;
    int windowsize = 30;
    bool loci = false;

    if (argc == 1)
        printSummary(argv);

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {"bed",  required_argument, 0, 'b'},
            {"invert",  no_argument, 0, 'v'},
	    {"intersect-vcf", required_argument, 0, 'i'},
	    {"union-vcf", required_argument, 0, 'u'},
            {"contained",  no_argument, 0, 'c'},
            {"overlapping", no_argument, 0, 'o'},
	    {"window-size", required_argument, 0, 'w'},
	    {"reference", required_argument, 0, 'r'},
	    {"loci", no_argument, 0, 'l'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvclob:i:u:w:r:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

	    case 'w':
		windowsize = atoi(optarg);
		break;

            case 'b':
                bedFileName = string(optarg);
                break;

            case 'i':
		intersecting = true;
                vcfFileName = string(optarg);
                break;

            case 'u':
		unioning = true;
                vcfFileName = string(optarg);
                break;

	    case 'r':
		fastaFileName = string(optarg);
		break;

            case 'v':
                invert = true;
                break;

            case 'c':
                contained = true;
                break;

            case 'o':
                overlapping = true;
                break;

	    case 'l':
	        loci = true;
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

    bool usingBED = false;
    if (!bedFileName.empty()) {
	usingBED = true;
    }
    
    BedReader bed;
    if (usingBED) {
	bed.open(bedFileName);
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

    if (usingBED) {
	variantFile.parseSamples = false;
    }

    VariantCallFile otherVariantFile;
    if (!vcfFileName.empty()) {
	otherVariantFile.open(vcfFileName);
	if (!otherVariantFile.is_open()) {
	    cerr << "could not open VCF file " << vcfFileName << endl;
	    exit(1);
	}
    }

    FastaReference reference;
    if (unioning || intersecting) {
	if (fastaFileName.empty()) {
	    cerr << "a reference is required for haplotype-based intersection and unioniong" << endl;
	    exit(1);
	}
	reference.open(fastaFileName);
    }

    if (!unioning && !intersecting) {
	variantFile.parseSamples = false; // faster, as when we are
					  // only bed-intersecting we
					  // can do position-only
					  // output and don't have to
					  // manipulate specific
					  // alleles
    }

    // read the VCF file for union or intersection into an interval tree
    // indexed using some proximity window

    map<string, IntervalTree<Variant*> > variantIntervals;
    map<string, list<Variant> > otherVariants;
    map<string, vector<Interval<Variant*> > > otherVariantIntervals;

    if (unioning || intersecting) {

	Variant ovar(otherVariantFile);
	while (otherVariantFile.getNextVariant(ovar)) {
	    long int left = ovar.position;
	    long int right = left + ovar.ref.size(); // this should be 1-past the end
	    otherVariants[ovar.sequenceName].push_back(ovar);
	    Variant* v = &otherVariants[ovar.sequenceName].back();
	    otherVariantIntervals[ovar.sequenceName].push_back(Interval<Variant*>(left, right, v));
	}
	
	for (map<string, vector<Interval<Variant*> > >::iterator j = otherVariantIntervals.begin(); j != otherVariantIntervals.end(); ++j) {
	    variantIntervals[j->first] = IntervalTree<Variant*>(j->second);
	}

    }

    set<Variant*> outputVariants;

    long unsigned int lastOutputPosition = 0;
    string lastSequenceName;

    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

	if (lastSequenceName.empty()) {
	    lastSequenceName = var.sequenceName;
	} else if (lastSequenceName != var.sequenceName) {
	    if (unioning) {
		vector<Interval<Variant*> > previousRecords;
		long int lastSeqLength = reference.sequenceLength(lastSequenceName);
		variantIntervals[lastSequenceName].findContained(lastOutputPosition, lastSeqLength, previousRecords);
		for (vector<Interval<Variant*> >::iterator r = previousRecords.begin(); r != previousRecords.end(); ++r) {
		    Variant* v = r->value;
		    if (outputVariants.find(v) == outputVariants.end()) {
			outputVariants.insert(v);
			cout << *v << endl;  // does this output everything in correct order?
		    }
		}
		lastSequenceName = var.sequenceName;
		lastOutputPosition = 0;
	    }
	}

	if (usingBED) {
	    BedTarget record(var.sequenceName, var.position, var.position + var.ref.size(), "");
	    vector<BedTarget*> overlaps = bed.targetsOverlapping(record);

	    if (!invert && !overlaps.empty()) {
		cout << variantFile.line << endl;
	    } else if (invert && overlaps.empty()) {
		cout << variantFile.line << endl;
	    }

	} else if (unioning || intersecting) {

	    // TODO check overlaps with union/intersection
	    // hmm... for unioning, you might need to step through the original VCF records
	    // but the idea is to exclude the haplotype-based duplicates

	    vector<Interval<Variant*> > results;

	    variantIntervals[var.sequenceName].findContained(var.position - windowsize, var.position + var.ref.size() + windowsize, results);

	    vector<Variant*> overlapping;

	    for (vector<Interval<Variant*> >::iterator r = results.begin(); r != results.end(); ++r) {
		overlapping.push_back(r->value);
	    }


	    if (unioning) {

		// unioning strategy

		// write out all the records from the last file
		// between the last one printed out and the first
		// one we're about to print out

		vector<Interval<Variant*> > previousRecords;

		variantIntervals[var.sequenceName].findOverlapping(lastOutputPosition, var.position - windowsize, previousRecords);

		map<long int, vector<Variant*> > variants;

		for (vector<Interval<Variant*> >::iterator r = previousRecords.begin(); r != previousRecords.end(); ++r) {
		    Variant* v = r->value;
		    if (outputVariants.find(v) == outputVariants.end()) {
			outputVariants.insert(v);
			variants[v->position].push_back(v);
		    }
		}

		for (map<long int, vector<Variant*> >::iterator v = variants.begin(); v != variants.end(); ++v) {
		    for (vector<Variant*>::iterator o = v->second.begin(); o != v->second.end(); ++o) {
			cout << **o << endl;
			lastOutputPosition = max(lastOutputPosition, (*o)->position);
		    }
		}

		// TODO find the duplicates for the other file
	    }


	    if (overlapping.empty()) {

		if (unioning || (intersecting && invert)) {
		    cout << var << endl;
		    lastOutputPosition = max(lastOutputPosition, var.position);
		}

	    } else {

		// get the min and max of the overlaps

		int haplotypeStart = var.position;
		int haplotypeEnd = var.position + var.ref.size();

		for (vector<Variant*>::iterator v = overlapping.begin(); v != overlapping.end(); ++v) {
		    haplotypeStart = min((*v)->position, (long unsigned int) haplotypeStart);
		    haplotypeEnd = max((*v)->position + (*v)->ref.size(), (long unsigned int) haplotypeEnd);
     		}

		// for everything overlapping and the current variant, construct the local haplotype within the bounds
		// if there is an exact match, the alllele in the current VCF does intersect

		string referenceHaplotype = reference.getSubSequence(var.sequenceName, haplotypeStart - 1, haplotypeEnd - haplotypeStart);
		map<string, vector<Variant*> > haplotypes;

		for (vector<Variant*>::iterator v = overlapping.begin(); v != overlapping.end(); ++v) {
		    Variant& variant = **v;
		    for (vector<string>::iterator a = variant.alt.begin(); a != variant.alt.end(); ++a) {
			string haplotype = referenceHaplotype;
			// get the relative start and end coordinates for the variant alternate allele
			int relativeStart = variant.position - haplotypeStart;
			haplotype.replace(relativeStart, variant.ref.size(), *a);
			haplotypes[haplotype].push_back(*v);
		    }
		}

		// determine the non-intersecting alts
		vector<string> altsToRemove;
		for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
		    string haplotype = referenceHaplotype;
		    int relativeStart = var.position - haplotypeStart;
		    haplotype.replace(relativeStart, var.ref.size(), *a);
		    map<string, vector<Variant*> >::iterator h = haplotypes.find(haplotype);
		    if ((intersecting && !invert && h == haplotypes.end())
			|| (intersecting && invert && h != haplotypes.end())
			|| (unioning && h != haplotypes.end())) {
			altsToRemove.push_back(*a);
		    }
		}

		// remove the non-overlapping (intersecting) or overlapping (unioning) alts
		if (intersecting && loci && altsToRemove.size() != var.alt.size()) {
		    // we have a match in loci mode, so we should output the whole loci, not just the matching sequence
		} else {
		    for (vector<string>::iterator a = altsToRemove.begin(); a != altsToRemove.end(); ++a) {
			var.removeAlt(*a);
		    }
		}

		if (unioning) {
		    // somehow sort the records and combine them?
		    map<long int, vector<Variant*> > variants;
		    for (vector<Variant*>::iterator o = overlapping.begin(); o != overlapping.end(); ++o) {
			if ((*o)->position <= var.position && // check ensures proper ordering of variants on output
			    outputVariants.find(*o) == outputVariants.end()) {
			    outputVariants.insert(*o);
			    variants[(*o)->position].push_back(*o);
			}
		    }
		    // add in the current variant, if it has alts left
		    if (!var.alt.empty()) {
			variants[var.position].push_back(&var);
		    }

		    for (map<long int, vector<Variant*> >::iterator v = variants.begin(); v != variants.end(); ++v) {
			for (vector<Variant*>::iterator o = v->second.begin(); o != v->second.end(); ++o) {
			    cout << **o << endl;
			    lastOutputPosition = max(lastOutputPosition, (*o)->position);
			}
		    }
		} else {
		    // if any alts remain, output the variant record
		    if (!var.alt.empty()) {
			cout << var << endl;
			lastOutputPosition = max(lastOutputPosition, var.position);
		    }
		}

	    }

	}

    }


    // if unioning, and any variants remain, output them
    if (unioning) {
	for (map<string, list<Variant> >::iterator chrom = otherVariants.find(lastSequenceName);
	     chrom != otherVariants.end();
	     ++chrom) {
	    for (list<Variant>::iterator v = chrom->second.begin(); v != chrom->second.end(); ++v) {
		Variant* variant = &*v;
		if (outputVariants.find(variant) == outputVariants.end()) {
		    outputVariants.insert(variant);
		    cout << *variant << endl;
		    // TODO guarantee sorting
		}
	    }
	}
    }

    exit(0);  // why?
    return 0;

}

