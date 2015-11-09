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
         << "    -b, --bed FILE            use intervals provided by this BED file" << endl
         << "    -R, --region REGION       use 1-based tabix-style region (e.g. chrZ:10-20), multiples allowed" << endl
         << "    -S, --start-only          don't use the reference length information in the record to determine" << endl
         << "                              overlap status, just use the start posiion" << endl
         << "    -v, --invert              invert the selection, printing only records which would" << endl
         << "                                not have been printed out" << endl
         << "    -i, --intersect-vcf FILE  use this VCF for set intersection generation" << endl
         << "    -u, --union-vcf FILE      use this VCF for set union generation" << endl
         << "    -w, --window-size N       compare records up to this many bp away (default 30)" << endl
         << "    -r, --reference FILE      FASTA reference file, required with -i and -u" << endl
         << "    -l, --loci                output whole loci when one alternate allele matches" << endl
         << "    -m, --ref-match           intersect on the basis of record REF string" << endl
         << "    -t, --tag TAG             attach TAG to each record's info field if it would intersect" << endl
         << "    -V, --tag-value VAL       use this value to indicate that the allele is passing" << endl
         << "                              '.' will be used otherwise.  default: 'PASS'" << endl
         << "    -M, --merge-from FROM-TAG" << endl
         << "    -T, --merge-to   TO-TAG   merge from FROM-TAG used in the -i file, setting TO-TAG" << endl
         << "                              in the current file." << endl
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
    bool startPositionOnly = false;
    int windowsize = 30;
    bool loci = false;
    bool refmatch = false;
    string tag;
    string tagValue = "PASS";
    string mergeFromTag;
    string mergeToTag;
    vector<BedTarget> regions;

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
                {"region",  required_argument, 0, 'R'},
                {"invert",  no_argument, 0, 'v'},
                {"intersect-vcf", required_argument, 0, 'i'},
                {"union-vcf", required_argument, 0, 'u'},
                {"contained",  no_argument, 0, 'c'},
                {"overlapping", no_argument, 0, 'o'},
                {"window-size", required_argument, 0, 'w'},
                {"reference", required_argument, 0, 'r'},
                {"loci", no_argument, 0, 'l'},
                {"ref-match", no_argument, 0, 'm'},
                {"tag", required_argument, 0, 't'},
                {"tag-value", required_argument, 0, 'V'},
                {"merge-from", required_argument, 0, 'M'},
                {"merge-to", required_argument, 0, 'T'},
                {"start-only", no_argument, 0, 'S'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hvcSlmob:i:u:w:r:t:V:M:T:R:",
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

	    case 'm':
	        refmatch = true;
	        break;

	    case 't':
	        tag = optarg;
            break;

        case 'R':
            regions.push_back(BedTarget(optarg));
            regions.back().left -= 1;
            break;

        case 'S':
            startPositionOnly = true;
            break;

	    case 'V':
            tagValue = optarg;
            break;

	    case 'M':
            mergeFromTag = optarg;
            break;

	    case 'T':
            mergeToTag = optarg;
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


    bool usingBED = false;
    if (!bedFileName.empty()) {
        usingBED = true;
    }
    
    if (usingBED || !regions.empty()) {
        variantFile.parseSamples = false;
    }

    // it runs much faster to do this first.  then downstream processes don't block!

    BedReader bed;
    if (usingBED) {
        bed.open(bedFileName);
    }
    if (!regions.empty()) {
        // add to the bed
        bed.addTargets(regions);
        usingBED = true;
    }

    VariantCallFile otherVariantFile;
    if (!vcfFileName.empty()) {
        if (vcfFileName == "-") {
            if (usingstdin) {
                cerr << "cannot open both VCF file streams from stdin" << endl;
                exit(1);
            } else {
                otherVariantFile.open(std::cin);
            }
        } else {
            otherVariantFile.open(vcfFileName);
        }
        if (!otherVariantFile.is_open()) {
            cerr << "could not open VCF file " << vcfFileName << endl;
            exit(1);
        }
    }


    if (!tag.empty()) {
        variantFile.addHeaderLine("##INFO=<ID="+ tag +",Number=A,Type=String,Description=\"" + tagValue + " if this allele intersects with one in " + vcfFileName  +  ", '.' if not.\">");
    }

    if (!mergeToTag.empty()) {
        if (mergeFromTag.empty()) {
            cerr << "must specify a tag to merge from" << endl;
            exit(1);
        }
        // get mergeFromTag type
        map<string, VariantFieldType>::iterator f = otherVariantFile.infoTypes.find(mergeFromTag);
        if (f == otherVariantFile.infoTypes.end()) {
            cerr << "vcfintersect: ERROR could not find " << mergeFromTag << " in header" << endl;
            exit(1);
        }
        VariantFieldType mergeFromType = f->second;
        stringstream s;
        s << mergeFromType;
        
        variantFile.addHeaderLine("##INFO=<ID="+ mergeToTag +",Number=A,Type=" + s.str() + ",Description=\"The value of " + mergeFromTag + " in " + vcfFileName  +  " '.' if the tag does not exist for the given allele in the other file, or if there is no corresponding allele.\">");
    }

    cout << variantFile.header << endl;


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

    long int lastOutputPosition = 0;
    string lastSequenceName;

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
                        cout << *v << endl; // Q: does this output everything in correct order?.... A: No.
                    }
                }
                lastSequenceName = var.sequenceName;
                lastOutputPosition = 0;
            }
        }

        if (usingBED) {
            vector<BedTarget*> overlaps;
            if (startPositionOnly) {
                // only intersect if start position (not end) is in target
                BedTarget record(var.sequenceName, var.position, var.position, "");
                overlaps = bed.targetsOverlapping(record);
            } else {
                // default behavior
                BedTarget record(var.sequenceName, var.position, var.position + var.ref.size() - 1, "");
                overlaps = bed.targetsOverlapping(record);
            }

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
                } else if (intersecting && (!tag.empty() || !mergeToTag.empty())) {
                    for (int i = 0; i < var.alt.size(); ++i) {
                        if (!tag.empty()) {
                            var.info[tag].push_back(".");
                        }
                        if (!mergeToTag.empty()) {
                            var.info[mergeToTag].push_back(".");
                        }
                    }
                    cout << var << endl;
                    lastOutputPosition = max(lastOutputPosition, var.position);
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

                Variant originalVar = var;

                // determine the non-intersecting alts
                vector<string> altsToRemove;
                vector<int> altIndexesToRemove;
                for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                    string haplotype = referenceHaplotype;
                    int relativeStart = var.position - haplotypeStart;
                    haplotype.replace(relativeStart, var.ref.size(), *a);
                    map<string, vector<pair<Variant*, int> > >::iterator h = haplotypes.find(haplotype);
                    if ((intersecting && !invert && h == haplotypes.end())
                        || (intersecting && invert && h != haplotypes.end())
                        || (unioning && h != haplotypes.end())) {
                        if (tag.empty() && mergeToTag.empty()) {
                            altsToRemove.push_back(*a);
                        } else {
                            if (!tag.empty()) {
                                var.info[tag].push_back(".");
                            }
                            if (!mergeToTag.empty()) {
                                var.info[mergeToTag].push_back(".");
                            }
                        }
                    } else {
                        if (!tag.empty()) {
                            var.info[tag].push_back(tagValue);
                        }
                        // NB: just take the first value for the mergeFromTag
                        if (!mergeToTag.empty()) {
                            Variant* v = h->second.front().first;
                            int index = h->second.front().second;
                            if (v->info.find(mergeFromTag) != v->info.end()) {
                                // now you have to find the exact allele...
                                string& otherValue = v->info[mergeFromTag].at(index);
                                var.info[mergeToTag].push_back(otherValue);
                            } else if (mergeFromTag == "QUAL") {
                                var.info[mergeToTag].push_back(convert(v->quality));
                            } else {
                                var.info[mergeToTag].push_back(".");
                            }
                        }
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
                        vector<Variant*>& vars = variants[var.position];
                        int numalts = 0;
                        for (vector<Variant*>::iterator v = vars.begin(); v != vars.end(); ++v) {
                            numalts += (*v)->alt.size();
                        }
                        if (numalts + var.alt.size() == originalVar.alt.size()) {
                            variants[var.position].clear();
                            variants[var.position].push_back(&originalVar);
                        } else {
                            variants[var.position].push_back(&var);
                        }
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

