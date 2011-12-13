#include "Variant.h"
#include "convert.h"
#include <set>

using namespace std;
using namespace vcf;


double convertStrDbl(const string& s) {
    double r;
    convert(s, r);
    return r;
}

int main(int argc, char** argv) {

    int maxAlleles = 2;

    VariantCallFile variantFile;

    if (argc > 1) {
        string filename = argv[1];
        if (filename == "--help" || filename == "-h") {
            cerr << "usage: " << argv[0] << " [file]" << endl
                 << endl
                 << "If multiple alleleic primitives (gaps or mismatches) are specified in" << endl
                 << "a single VCF record, split the record into multiple lines, but drop all" << endl
                 << "INFO fields.  Does not handle genotypes (yet)." << endl;
            exit(1);
        }
        variantFile.open(filename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        return 1;
    }

    cout << variantFile.header;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

        //cout << var << endl;

        // for each parsedalternate, get the position
        // build a new vcf record for that position
        // unless we are already at the position !
        // take everything which is unique to that allele (records) and append it to the new record
        bool includePreviousBaseForIndels = true;
        map<string, vector<VariantAllele> > varAlleles = var.parsedAlternates(includePreviousBaseForIndels);
        set<VariantAllele> alleles;
        for (map<string, vector<VariantAllele> >::iterator a = varAlleles.begin(); a != varAlleles.end(); ++a) {
            for (vector<VariantAllele>::iterator va = a->second.begin(); va != a->second.end(); ++va) {
                alleles.insert(*va);
            }
        }

        //map<long unsigned int, Variant> variants;
        vector<Variant> variants;
        for (set<VariantAllele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            if (a->ref == a->alt) {
                // ref allele
                continue;
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

        //for (map<long unsigned int, Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
        for (vector<Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
            cout << *v << endl;
        }
    }

    return 0;

}

