#include "Variant.h"
#include <getopt.h>

using namespace std;
using namespace vcf;

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [vcf file] [vcf file] ..." << endl
         << endl
         << "Combines VCF files positionally, combining samples when sites and alleles are identical." << endl
         << "Any number of VCF files may be combined.  The INFO field and other columns are taken from" << endl
         << "one of the files which are combined when records in multiple files match.  Alleles must" << endl
         << "have identical ordering to be combined into one record.  If they do not, multiple records" << endl
         << "will be emitted." << endl;
    exit(1);
}


int main(int argc, char** argv) {

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "h",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

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

    vector<string> sampleNames;
    string randomHeader;
    VariantCallFile* vcf;

    list<string> chromOrder;

    map<string, // chrom
        map<long int, // pos
            map<vector<string>, // alts
                map<VariantCallFile*, Variant*> > > > variantFiles;
    for (int i = optind; i != argc; ++i) {
        string inputFilename = argv[i];
        vcf = new VariantCallFile;
        vcf->open(inputFilename);
        Variant* var = new Variant(*vcf);
        vcf->getNextVariant(*var);
        if (variantFiles.find(var->sequenceName) == variantFiles.end()) {
            chromOrder.push_back(var->sequenceName);
        }
        variantFiles[var->sequenceName][var->position][var->alt][vcf] = var;
        sampleNames.insert(sampleNames.end(), vcf->sampleNames.begin(), vcf->sampleNames.end());
    }

    sort(sampleNames.begin(), sampleNames.end());
    sampleNames.erase(unique(sampleNames.begin(), sampleNames.end()), sampleNames.end());

    VariantCallFile outputCallFile;
    string header = vcf->headerWithSampleNames(sampleNames);
    outputCallFile.openForOutput(header);

    cout << outputCallFile.header << endl;

    while (!variantFiles.empty()) {
        // get lowest variant(s)
        // if they have identical alts and position, combine
        // otherwise just output, but with the same sample names
        map<long int, map<vector<string>, map<VariantCallFile*, Variant*> > >& chrom = variantFiles[chromOrder.front()];
        if (chrom.empty()) {
            variantFiles.erase(chromOrder.front());
            chromOrder.pop_front();
            continue;
        }
        map<vector<string>, map<VariantCallFile*, Variant*> >& pos = chrom.begin()->second;
        map<vector<string>, map<VariantCallFile*, Variant*> >::iterator s = pos.begin();
        for ( ; s != pos.end(); ++s) {
            Variant variant(outputCallFile);
            map<VariantCallFile*, Variant*>& vars = s->second;
            map<VariantCallFile*, Variant*>::iterator v = vars.begin();
            for ( ; v != vars.end(); ++v) {
                VariantCallFile* vcf = v->first;
                Variant* var = v->second;
                if (variant.info.empty()) {
                    variant.sequenceName = var->sequenceName;
                    variant.position = var->position;
                    variant.id = var->id;
                    variant.ref = var->ref;
                    variant.alt = var->alt;
                    variant.quality = var->quality;
                    variant.info = var->info;
                    variant.format = var->format;
                }
                // add samples to output variant
                for (Samples::iterator sample = var->samples.begin(); sample != var->samples.end(); ++sample) {
                    variant.samples[sample->first] = sample->second;
                }
                if (vcf->getNextVariant(*var)) {
                    if (variantFiles.find(var->sequenceName) == variantFiles.end()) {
                        chromOrder.push_back(var->sequenceName);
                    }
                    variantFiles[var->sequenceName][var->position][var->alt][vcf] = var;
                }
            }
            if (!variant.info.empty())
                cout << variant << endl;
        }
        // XXXXXXX this will cause problems if the files have more than one record per position!!!
        chrom.erase(chrom.begin());
    }

    return 0;

}

