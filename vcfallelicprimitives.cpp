#include "Variant.h"
#include "convert.h"
#include <set>
#include <getopt.h>

using namespace std;
using namespace vcf;


double convertStrDbl(const string& s) {
    double r;
    convert(s, r);
    return r;
}

void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [file]" << endl
	 << endl
	 << "If multiple alleleic primitives (gaps or mismatches) are specified in" << endl
	 << "a single VCF record, split the record into multiple lines, but drop all" << endl
	 << "INFO fields.  Does not handle genotypes (yet).  MNPs are split into multiple SNPs." << endl;
    exit(0);
}

int main(int argc, char** argv) {

    bool includePreviousBaseForIndels = true;
    bool useMNPs = false;

    VariantCallFile variantFile;

    int c;
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"help", no_argument, 0, 'h'},
	    {"use-mnps", no_argument, 0, 'm'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hm",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

	    case 'm':
		useMNPs = true;
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

    if (optind < argc) {
        string filename = argv[1];
        variantFile.open(filename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    variantFile.addHeaderLine("##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">");
    variantFile.addHeaderLine("##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">");
    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

        //cout << var << endl;

        // for each parsedalternate, get the position
        // build a new vcf record for that position
        // unless we are already at the position !
        // take everything which is unique to that allele (records) and append it to the new record

        map<string, vector<VariantAllele> > varAlleles = var.parsedAlternates(includePreviousBaseForIndels, useMNPs);
        set<VariantAllele> alleles;

        for (map<string, vector<VariantAllele> >::iterator a = varAlleles.begin(); a != varAlleles.end(); ++a) {
            for (vector<VariantAllele>::iterator va = a->second.begin(); va != a->second.end(); ++va) {
                alleles.insert(*va);
            }
        }

	map<VariantAllele, double> alleleFrequencies;
	map<VariantAllele, int> alleleCounts;

	bool hasAf = false;
	if (var.info.find("AF") != var.info.end()) {
	    hasAf = true;
	    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
		vector<VariantAllele>& vars = varAlleles[*a];
		for (vector<VariantAllele>::iterator va = vars.begin(); va != vars.end(); ++va) {
		    double freq;
		    convert(var.info["AF"].at(var.altAlleleIndexes[*a]), freq);
		    alleleFrequencies[*va] += freq;
		}
	    }
	}

	bool hasAc = false;
	if (var.info.find("AC") != var.info.end()) {
	    hasAc = true;
	    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
		vector<VariantAllele>& vars = varAlleles[*a];
		for (vector<VariantAllele>::iterator va = vars.begin(); va != vars.end(); ++va) {
		    int freq;
		    convert(var.info["AC"].at(var.altAlleleIndexes[*a]), freq);
		    alleleCounts[*va] += freq;
		}
	    }
	}

        //map<long unsigned int, Variant> variants;
        vector<Variant> variants;
        for (set<VariantAllele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            if (a->ref == a->alt) {
                // ref allele
                continue;
            }
	    string type;
	    int len = 0;
	    if (a->ref.at(0) == a->alt.at(0)) { // well-behaved indels
		if (a->ref.size() > a->alt.size()) {
		    type = "del";
		    len = a->ref.size() - a->alt.size();
		} else if (a->ref.size() < a->alt.size()) {
		    len = a->alt.size() - a->ref.size();
		    type = "ins";
		}
	    } else {
		if (a->ref.size() == a->alt.size()) {
		    len = a->ref.size();
		    if (a->ref.size() == 1) {
			type = "snp";
		    } else {
			type = "mnp";
		    }
		} else {
		    len = abs((int) a->ref.size() - (int) a->alt.size());
		    type = "complex";
		}
	    }
            //cout << a->ref << "/" << a->alt << endl;
            /*
            if (variants.find(a->position) == variants.end()) {
                Variant newvar(variantFile);
                newvar.quality = var.quality;
                newvar.filter = ".";
                variants.insert(make_pair(a->position, newvar));
            }
            Variant& v = variants[a->position];
            */
            Variant newvar(variantFile);
            newvar.quality = var.quality;
            newvar.filter = ".";
            newvar.id = ".";
	    newvar.format = var.format;
	    newvar.info["TYPE"].push_back(type);
	    newvar.info["LEN"].push_back(convert(len));
	    if (hasAf) {
		newvar.info["AF"].push_back(convert(alleleFrequencies[*a]));
	    }
	    if (hasAc) {
		newvar.info["AC"].push_back(convert(alleleCounts[*a]));
	    }
            variants.push_back(newvar);
            Variant& v = variants.back();
            v.sequenceName = var.sequenceName;
            v.position = a->position;
            if (v.ref.size() < a->ref.size()) {
                for (vector<string>::iterator va = v.alt.begin(); va != v.alt.end(); ++va) {
                    *va += a->ref.substr(v.ref.size());
                }
                v.ref = a->ref;
            }
            v.alt.push_back(a->alt);
        }

	// TODO genotypes!

        //for (map<long unsigned int, Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
        for (vector<Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
            cout << *v << endl;
        }
    }

    return 0;

}

