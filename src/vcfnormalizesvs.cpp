#include "Variant.h"
#include <getopt.h>
#include "Fasta.h"

using namespace std;
using namespace vcflib;

void print_help(char** argv){
    cerr << "vcfnormalizesvs: Convert a VCF with representative alleles to one with canonical, sequence-based ref/alt alleles." << endl
        << "usage: " << argv[0] << " [OPTIONS] var.vcf" << endl
        << "Options: " << endl
        << "   -r / --reference <ref.fa>   FASTA-format reference genome from which to pull SV sequences." << endl
        << "   -i / --inserions <ins.fa>   FASTA-format insertion sequences, with IDs matching the ALT allele tags in the vcf" << endl
        << endl;
}

int main(int argc, char** argv) {

    
    string ref_file = "";
    vector<string> insertion_files;
    int max_interval = -1;

    int c = 0;
    while (true) {
        static struct option long_options[] =
            {
                {"insertions", no_argument, 0, 'i'},
                {"help", no_argument, 0, 'h'},
                {"reference", required_argument, 0, 'r'},
                {0, 0, 0, 0}
            };
        int option_index = 0;

        c = getopt_long (argc, argv, "r:i:h",
                         long_options, &option_index);
        if (c == -1)
            break;
        /* Detect the end of the options. */
        switch(c){
        case 'r':
            ref_file = optarg;
            break;
        case 'i':
            insertion_files.push_back(optarg);
            break;
        case 'h':
        case '?':
            print_help(argv);
            exit(1);
        default:
            print_help(argv);
            abort();
        }
    }

    if (argc < 2){
        print_help(argv);
        exit(1);
    }

 

    VariantCallFile variantFile;
    string filename = argv[argc - 1];
    variantFile.open(filename);
    if (!variantFile.is_open()) {
        return 1;
    }

    vector<FastaReference*> insertions;
    if (!insertion_files.empty()){
        for (auto x : insertion_files){
            FastaReference* ins = new FastaReference();
            insertions.push_back(ins);
            ins->open(x);
        }
    }

    FastaReference ref;
    if(!ref_file.empty()){
        ref.open(ref_file);
    }


    cout << variantFile.header << endl;

    Variant var;
    while (variantFile.getNextVariant(var)) {
        bool valid = var.canonicalize_sv(ref, insertions, max_interval);
        if (!valid){
            cerr << "Variant could not be normalized" << var << endl;
        }
        cout << var << endl;
    }

    return 0;

}

