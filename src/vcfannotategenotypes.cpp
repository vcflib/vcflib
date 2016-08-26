#include "Variant.h"
#include "split.h"
#include <string>
#include <list>
#include <iostream>

using namespace std;
using namespace vcflib;

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
           
            if (otherGenotype.find("|") != string::npos) {
              vector<int> gtB = decomposePhasedGenotype(otherGenotype);
              vector<int> gtnew;
              gtnew.reserve(gtB.size());

              for (vector<int>::iterator g = gtB.begin(); g != gtB.end(); ++g) {
                map<int, int>::iterator f = varBconvertToVarA.find(*g);
                if (f != varBconvertToVarA.end()) {
                  gtnew.push_back(*g);
                } else {
                  gtnew.push_back(NULL_ALLELE);
                }
              }
              sample[annotationTag].push_back(phasedGenotypeToString(gtnew));
            } else {
              map<int, int> gtB = decomposeGenotype(otherGenotype);
              map<int, int> gtnew;
              for (map<int, int>::iterator g = gtB.begin(); g != gtB.end(); ++g) {
                map<int, int>::iterator f = varBconvertToVarA.find(g->first);
                if (f != varBconvertToVarA.end()) {
                  gtnew[f->second] += g->second;
                } else {
                  gtnew[NULL_ALLELE] += g->second;
                }
              }
              sample[annotationTag].push_back(genotypeToString(gtnew));
            }
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
    line = "##FORMAT=<ID=" + annotag + ",Number=1,Type=String,Description=\"Genotype from "
        + annotag + ".\">";
    variantFileA.addHeaderLine(line);

    cout << variantFileA.header << endl;

    do {

        // this is broken.  to do it right, it'll be necessary to get reference ids from the fasta reference used to make the alignments...
		// if B is NOT done, and is less than A, read new B.
        if (!variantFileB.done()
            && (varB.sequenceName != varA.sequenceName
                || (varB.sequenceName == varA.sequenceName && varB.position < varA.position)
				|| variantFileA.done())
            ) {
            variantFileB.getNextVariant(varB);
        }

		// if A is not done- and A is less than B, read A.  
		// should also read if variant B is done. 
        if (!variantFileA.done()
            && (varA.sequenceName != varB.sequenceName
                || (varA.sequenceName == varB.sequenceName && varA.position < varB.position)
				|| variantFileB.done())
            ) {
            annotateWithBlankGenotypes(varA, annotag);
            cout << varA << endl;
            variantFileA.getNextVariant(varA);
        }

        vector<Variant> varsA;
        vector<Variant> varsB;

        bool hasMultipleAlts = false;

        long int thisPosition = 0;
        string thisSequenceName;
        if (varA.position == varB.position
            && varA.sequenceName == varB.sequenceName) {
            thisPosition = varA.position;
            thisSequenceName = varA.sequenceName;
        }
        while (!variantFileA.done()
               && !variantFileB.done()
               && thisPosition == varA.position
               && thisSequenceName == varA.sequenceName
               && varA.sequenceName == varB.sequenceName
               && varA.position == varB.position) {
            // accumulate all the alts at the current position
            varsA.push_back(varA);
            varsB.push_back(varB);
            if (varA.alt.size() > 1 || varB.alt.size() > 1)
                hasMultipleAlts = true;
            variantFileA.getNextVariant(varA);
            variantFileB.getNextVariant(varB);
        }

        // multiple lines per position
        if (!hasMultipleAlts && (varsA.size() > 1 || varsB.size() > 1)) {

            map<pair<string, string>, Variant> varsAParsed;
            map<pair<string, string>, Variant> varsBParsed;	
            for (vector<Variant>::iterator v = varsA.begin(); v != varsA.end(); ++v) {
                varsAParsed[make_pair(v->ref, v->alt.front())] = *v;
            }
            for (vector<Variant>::iterator v = varsB.begin(); v != varsB.end(); ++v) {
                varsBParsed[make_pair(v->ref, v->alt.front())] = *v;
            }
	    
            for (map<pair<string, string>, Variant>::iterator vs = varsAParsed.begin(); vs != varsAParsed.end(); ++vs) {
                Variant& varA = vs->second;
                annotateWithBlankGenotypes(varA, annotag);
                if (varsBParsed.find(make_pair(varA.ref, varA.alt.front())) != varsBParsed.end()) {
                    Variant& varB = varsBParsed[make_pair(varA.ref, varA.alt.front())]; // TODO cleanup
                    annotateWithGenotypes(varA, varB, annotag);
                    varA.infoFlags[annotag + ".has_variant"] = true;
                }
                cout << varA << endl;
            }

        } else if (!varsA.empty() && !varsB.empty()) { // one line per multi-allelic
            Variant& varA = varsA.front();
            annotateWithBlankGenotypes(varA, annotag);
            Variant& varB = varsB.front();
            annotateWithGenotypes(varA, varB, annotag);
            // XXX TODO, and also allow for records with multiple alts
            // XXX assume that if the other file has a corresponding record, some kind of variation was detected at the same site
            varA.infoFlags[annotag + ".has_variant"] = true;
            cout << varA << endl;
        } else {
            for (vector<Variant>::iterator v = varsA.begin(); v != varsA.end(); ++v) {
                Variant& varA = *v;
                annotateWithBlankGenotypes(varA, annotag);
                cout << varA << endl;
            }
        }

    } while (!variantFileA.done() || !variantFileB.done());

    return 0;

}

