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
         << "    -r, --reference FILE         FASTA reference file." << endl
         << endl
         << "Modify the VCF file so that homozygous regions are included as REF/. calls." << endl
         << "For each ref and alt allele, assign an index.  These steps are sufficient to" << endl
         << "enable use of the VCF as a DAG (specifically a partially-ordered graph)." << endl;
    exit(0);
}

int main(int argc, char** argv) {

    string vcfFileName;
    string fastaFileName;

    bool adjustVcf = false;

    if (argc == 1)
        printSummary(argv);

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"reference", required_argument, 0, 'r'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hr:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

	    case 'r':
            fastaFileName = string(optarg);
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
        cerr << "a reference is required" << endl;
        exit(1);
    } else {
        reference.open(fastaFileName);
    }
    
    string idname = "id";
    long int uid = 0;

    variantFile.addHeaderLine("##INFO=<ID="+idname+".alt,Number=A,Type=Integer,Description=\"Unique numerical identifier of alt allele.\">");
    variantFile.addHeaderLine("##INFO=<ID="+idname+".ref,Number=1,Type=Integer,Description=\"Unique numerical identifier of ref allele.\">");
    cout << variantFile.header << endl;

    long int last_end = 1;
    string sequenceName;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

        if (sequenceName.empty()) {
            sequenceName = var.sequenceName;
        } else if (sequenceName != var.sequenceName) {
            // emit last record from previous chrom
            // these should be refactored.....
            Variant refvar(variantFile);
            if (var.position - last_end > 0) {
                refvar.ref = reference.getSubSequence(sequenceName, last_end - 1, var.position - last_end);
                refvar.quality = 0;
                refvar.position = last_end;
                refvar.sequenceName = sequenceName;
                refvar.info[idname+".ref"].push_back(convert(uid++));
                cout << refvar << endl;
            }
            last_end = 1;
            sequenceName = var.sequenceName;
        }

        // generate the last reference record if we have sequence between variants
        if (var.position - last_end > 0) {
            Variant refvar(variantFile);
            refvar.quality = 0;
            refvar.position = last_end;
            refvar.sequenceName = sequenceName;
            refvar.ref = reference.getSubSequence(sequenceName, last_end - 1, var.position - last_end);
            refvar.info[idname+".ref"].push_back(convert(uid++));
            cout << refvar << endl;
        }

        // now manipulate this record
        vector<string>& refidx = var.info[idname+".ref"];
        refidx.clear(); refidx.push_back(convert(uid++));

        vector<string>& idxs = var.info[idname+".alt"];
        idxs.clear();
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            idxs.push_back(convert(uid++));
        }
        cout << var << endl;

        last_end = var.position + var.ref.size();

    }

    if (reference.sequenceLength(sequenceName) - last_end > 0) {
        Variant refvar(variantFile);
        refvar.quality = 0;
        refvar.position = last_end;
        refvar.sequenceName = sequenceName;
        refvar.ref = reference.getSubSequence(sequenceName, last_end,
                                              reference.sequenceLength(sequenceName) - last_end);
        refvar.info[idname+".ref"].push_back(convert(uid++));
        cout << refvar << endl;
    }

    return 0;

}

