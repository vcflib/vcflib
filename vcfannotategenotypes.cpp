#include "Variant.h"
#include "split.h"
#include <string>
#include <iostream>

using namespace std;
using namespace vcf;

void annotateWithBlankGenotypes(Variant& var, string& annotationTag) {

    var.addFormatField(annotationTag);

    map<string, map<string, vector<string> > >::iterator s     = var.samples.begin(); 
    map<string, map<string, vector<string> > >::iterator sEnd  = var.samples.end();

    for (; s != sEnd; ++s) {
        map<string, vector<string> >& sample = s->second;
        sample[annotationTag].clear(); // means "no genotype" genotype
        sample[annotationTag].push_back("./."); // means "no genotype" genotype
    }
}

void annotateWithGenotypes(Variant& varA, Variant& varB, string& annotationTag) {

    varA.addFormatField(annotationTag);

    map<string, map<string, vector<string> > >::iterator s     = varA.samples.begin(); 
    map<string, map<string, vector<string> > >::iterator sEnd  = varA.samples.end();

    map<string, int> varAAlleleInts;
    int i = 1;
    for (vector<string>::iterator a = varA.alt.begin(); a != varA.alt.end(); ++a, ++i) {
        varAAlleleInts[*a] = i;
    }

    map<int, int> varBconvertToVarA; // maps alleles in the second file to allele numbers for the first
    varBconvertToVarA[0] = 0; // reference == reference!
    i = 1;
    for (vector<string>::iterator a = varB.alt.begin(); a != varB.alt.end(); ++a, ++i) {
        map<string, int>::iterator ita = varAAlleleInts.find(*a);
        if (ita != varAAlleleInts.end()) {
            varBconvertToVarA[i] = ita->second;
        }
    }

    for (; s != sEnd; ++s) {
        map<string, vector<string> >& sample = s->second;
        const string& name = s->first;
        map<string, map<string, vector<string> > >::iterator o = varB.samples.find(name);
        sample[annotationTag].clear();
        if (o == varB.samples.end()) {
            sample[annotationTag].push_back("./."); // means "no genotype"
        } else {
            map<string, vector<string> >& other = o->second;
            string& otherGenotype = other["GT"].front();
            // XXX this must compare the genotypes in the two files
            map<int, int> gtB = decomposeGenotype(otherGenotype);
            map<int, int> gtnew;
            for (map<int, int>::iterator g = gtB.begin(); g != gtB.end(); ++g) {
                map<int, int>::iterator f = varBconvertToVarA.find(g->first);
                if (f != varBconvertToVarA.end()) {
                    gtnew[f->second] += g->second;
                } else {
                    gtnew[-1] += g->second;
                }
            }
            sample[annotationTag].push_back(genotypeToString(gtnew));
        }
    }

}

int main(int argc, char** argv) {

    if (argc != 4) {
        cerr << "usage: " << argv[0] << " <annotation-tag> <vcf file> <vcf file>" << endl
             << "annotates genotypes in the first file with genotypes in the second" << endl
             << "adding the genotype as another flag to each sample filed in the first file." << endl
             << "annotation-tag is the name of the sample flag which is added to store the annotation." << endl
             << "also adds a 'has_variant' flag for sites where the second file has a variant." << endl;
        return 1;
    }

    string annotag = argv[1];
    string filenameA = argv[2];
    string filenameB = argv[3];

    if (filenameA == filenameB) {
        cerr << "it won't help to annotate samples with their own genotypes!" << endl;
        return 1;
    }

    VariantCallFile variantFileA;
    if (filenameA == "-") {
        variantFileA.open(std::cin);
    } else {
        variantFileA.open(filenameA);
    }

    VariantCallFile variantFileB;
    if (filenameB == "-") {
        variantFileB.open(std::cin);
    } else {
        variantFileB.open(filenameB);
    }

    if (!variantFileA.is_open() || !variantFileB.is_open()) {
        return 1;
    }

    Variant varA(variantFileA);
    Variant varB(variantFileB);

    // while the first file doesn't match the second positionally,
    // step forward, annotating each genotype record with an empty genotype
    // when the two match, iterate through the genotypes from the first file
    // and get the genotypes reported in the second file
    
    variantFileA.getNextVariant(varA);
    variantFileB.getNextVariant(varB);

    string line = "##INFO=<ID=" + annotag + ".has_variant,Number=0,Type=Flag,Description=\"True if "
        + annotag + " has a called alternate among samples under comparison.\">";
    variantFileA.addHeaderLine(line);

    cout << variantFileA.header << endl;

    if (variantFileA.done()) {
        cout << "is done" << endl;
    }

    if (variantFileB.done()) {
        cout << "is done" << endl;
    }


    do {

        while (!variantFileB.done() && varB.sequenceName < varA.sequenceName
                || (varB.sequenceName == varA.sequenceName && varB.position < varA.position)
                && !variantFileB.done()) {
            variantFileB.getNextVariant(varB);
        }

        while (!variantFileA.done()
                && varA.sequenceName < varB.sequenceName
                || (varA.sequenceName == varB.sequenceName && varA.position < varB.position)
                && !variantFileA.done()) {
            annotateWithBlankGenotypes(varA, annotag);
            cout << varA << endl;
            variantFileA.getNextVariant(varA);
        }

        while (!variantFileB.done() && varB.sequenceName < varA.sequenceName
                || (varB.sequenceName == varA.sequenceName && varB.position < varA.position)
                && !variantFileB.done()) {
            variantFileB.getNextVariant(varB);
        }

        while (!variantFileA.done() && varA.sequenceName == varB.sequenceName && varA.position == varB.position) {
            annotateWithGenotypes(varA, varB, annotag);
            // XXX assume that if the other file has a corresponding record, some kind of variation was detected at the same site
            varA.infoFlags[annotag + ".has_variant"] = true;
            cout << varA << endl;
            variantFileA.getNextVariant(varA);
        }
        
    } while (!variantFileA.done() && !variantFileB.done());

    return 0;

}

