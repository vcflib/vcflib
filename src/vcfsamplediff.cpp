#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>
#include <getopt.h>

using namespace std;
using namespace vcflib;

bool samplesDiffer(vector<string>& samples, Variant& var) {

    string genotype;

    for (vector<string>::iterator s = samples.begin(); s != samples.end(); ++s) {
        string& sampleName = *s;
        map<string, map<string, vector<string> > >::iterator f = var.samples.find(sampleName);
        if (f != var.samples.end()) {
            map<string, vector<string> >& sample = f->second;
            map<string, vector<string> >::iterator gt = sample.find("GT");
            if (gt != sample.end()) {
                string& thisGenotype = gt->second.front();
                if (genotype.empty()) {
                    genotype = thisGenotype;
                } else {
                    if (genotype != thisGenotype) {
                        return true;
                    }
                }
            }
        }
    }

    return false;

}


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <tag> <sample> <sample> [ <sample> ... ] <vcf file>" << endl
         << "Tags each record where the listed sample genotypes differ with <tag>." << endl
         << "The first sample is assumed to be germline, the second somatic." << endl
         << "Each record is tagged with <tag>={germline,somatic,loh} to specify the type of" << endl
         << "variant given the genotype difference between the two samples." << endl
         << endl
         << "options:" << endl
         << "    -s --strict     Require that no observations in the germline support the somatic alternate." << endl
         << endl;
}


int main(int argc, char** argv) {

    bool strict = false;
    int c;

    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"strict",  no_argument, 0, 's'},
                //{"length",  no_argument, &printLength, true},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hs",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
            printf ("\n");
            break;

        case 's':
            strict = true;
            break;
 
        case 'h':
            printSummary(argv);
            exit(0);
            break;

        case '?':
            /* getopt_long already printed an error message. */
            printSummary(argv);
            exit(1);
            break;
 
        default:
            abort ();
        }
    }

    if(argc - optind < 4) {
        printSummary(argv);
        exit(0);
    }

    string tag = argv[optind];

    vector<string> samples;
    for (int i = optind+1; i < argc - 1; ++i) {
        samples.push_back(argv[i]);
    }

    string filename = argv[argc-1];

    VariantCallFile variantFile;
    if (filename == "-") {
        variantFile.open(std::cin);
    } else {
        variantFile.open(filename);
    }

    if (!variantFile.is_open()) {
        cerr << "could not open " << filename << endl;
        return 1;
    }

    assert(samples.size() == 2);

    Variant var(variantFile);

    // TODO check if AC is present
    // ensure that AC is listed as an info field
    string line = "##INFO=<ID=" + tag + ",Number=1,Type=String,Description=\"Samples";
    for (vector<string>::iterator s = samples.begin(); s != samples.end(); ++s) {
        line += " " + *s;
    }
    line += " have different genotypes\">";
    variantFile.addHeaderLine(line);

    variantFile.addHeaderLine("##INFO=<ID=SSC,Number=1,Type=Float,Description=\"Somatic variant score (phred-scaled probability that the somatic variant call is correct).\">");

    // write the new header
    cout << variantFile.header << endl;
 
    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        if (var.samples.find(samples.front()) != var.samples.end()
            && var.samples.find(samples.back()) != var.samples.end()) {
            map<string, vector<string> >& germline = var.samples[samples.front()];
            map<string, vector<string> >& somatic = var.samples[samples.back()];
            map<int, int> gtGermline = decomposeGenotype(germline["GT"].front());
            map<int, int> gtSomatic  = decomposeGenotype(somatic["GT"].front());
            int germlineAltCount = 0;
            if (germline.find("AO") != germline.end()) {
                convert(germline["AO"].front(), germlineAltCount);
            }
            var.info[tag].clear(); // remove previous
            if (gtGermline == gtSomatic) {
                var.info[tag].push_back("germline");
            } else {
                //if (isHet(gtGermline) && isHom(gtSomatic)) {
                //    var.info[tag].push_back("loh");
                if (isHet(gtGermline) && isHomNonRef(gtSomatic) ||
                    isHomRef(gtGermline) && (isHet(gtSomatic) || isHomNonRef(gtSomatic))) {
                    if (!strict || strict && germlineAltCount == 0) {
                        var.info[tag].push_back("somatic");
                    }
                } else if (isHom(gtGermline) && isHet(gtSomatic)) {
                    if (var.alt.size() == 1) {
                        var.info[tag].push_back("reversion");
                    } else {
                        var.info[tag].push_back("somatic");
                    }
                }
            }
            if (germline.find("GQ") != germline.end() && somatic.find("GQ") != somatic.end()) {
                double germlineGQ;
                convert(germline["GQ"].front(), germlineGQ);
                double somaticGQ;
                convert(somatic["GQ"].front(), somaticGQ);
                double somaticScore = min(var.quality, min(germlineGQ, somaticGQ));
                var.info["SSC"].clear();
                var.info["SSC"].push_back(convert(somaticScore));
            }
        }
        cout << var << endl;
    }

    return 0;

}

