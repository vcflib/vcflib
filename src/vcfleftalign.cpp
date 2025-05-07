/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include "join.h"
#include "LeftAlign.hpp"

#include <Fasta.h>
#include <vector>
#include <getopt.h>
#include <cmath>
#include <SmithWatermanGotoh.h>


using namespace std;
using namespace vcflib;

void getAlignment(Variant& var, FastaReference& reference, string& ref, vector<AltAlignment>& alignments, int window) {

    // default alignment params
    float matchScore = 10.0f;
    float mismatchScore = -9.0f;
    float gapOpenPenalty = 25.0f;
    float gapExtendPenalty = 3.33f;

    // establish reference sequence
    string pad = string(window/2, 'Z');
    string leftFlank = reference.getSubSequence(var.sequenceName, var.zeroBasedPosition() - window/2, window/2);
    string rightFlank = reference.getSubSequence(var.sequenceName, var.zeroBasedPosition() + var.ref.size(), window/2);
    ref = pad + leftFlank + var.ref + rightFlank + pad;

    // and iterate through the alternates, generating alignments
    for (const auto& a : var.alt) {
        string alt = pad + leftFlank + a + rightFlank + pad;
        CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
        unsigned int referencePos;
        string cigar;
        sw.Align(referencePos, cigar, ref, alt);
        alignments.push_back(AltAlignment(referencePos, alt, cigar));
    }
}

void printSummary(char** argv) {
      cerr << R"(

Left-align indels and complex variants in the input using a pairwise
ref/alt alignment followed by a heuristic, iterative left realignment
process that shifts indel representations to their absolute leftmost
(5') extent.

This is the same procedure used in the internal left alignment in
freebayes, and can be used when preparing VCF files for input to
freebayes to decrease positional representation differences between
the input alleles and left-realigned alignments.

usage: vcfleftalign [options] [file]

options:

        -r, --reference FILE  Use this reference as a basis for realignment.
        -w, --window N        Use a window of this many bp when left aligning (150).

Left-aligns variants in the specified input file or stdin.  Window
size is determined dynamically according to the entropy of the regions
flanking the indel.  These must have entropy > 1 bit/bp, or be shorter
than ~5kb.

)";
    cerr << endl << "Type: transformation" << endl << endl;
    exit(0);
}

int main(int argc, char** argv) {

    int window = 150;
    VariantCallFile variantFile;
    string fastaFileName;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"help", no_argument, 0, 'h'},
                {"reference", required_argument, 0, 'r'},
                {"window", required_argument, 0, 'w'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hw:r:",
                         long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

	    case 'r':
            fastaFileName = optarg;
            break;

	    case 'w':
            window = atoi(optarg);
            break;

        case '?':
            printSummary(argv);
            exit(1);
            break;

        case 'h':
            printSummary(argv);
            break;

        default:
            abort ();
        }
    }

    if (optind < argc) {
        string filename = argv[optind];
        variantFile.open(filename);
    } else {
        variantFile.open(std::cin);
    }

    if (!variantFile.is_open()) {
        cerr << "could not open VCF file" << endl;
        exit(1);
    }

    FastaReference fastaReference;
    if (fastaFileName.empty()) {
        cerr << "a reference is required" << endl;
        exit(1);
    } else {
        fastaReference.open(fastaFileName);
    }

    /*
    variantFile.addHeaderLine("##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">");
    variantFile.addHeaderLine("##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">");
    if (!parseFlag.empty()) {
        variantFile.addHeaderLine("##INFO=<ID="+parseFlag+",Number=0,Type=Flag,Description=\"The allele was parsed using vcfallelicprimitives.\">");
    }
    */
    cout << variantFile.header << endl;

    Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {

        // if there is no indel, there is nothing to realign
        bool hasIndel = false;
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            if (a->size() != var.ref.size()) {
                hasIndel = true;
                break;
            }
        }
        if (!hasIndel) {
            cout << var << endl;
            continue;
        }

        vector<AltAlignment> alignments;
        string ref;

        // determine window size to prevent mismapping with SW algorithm
        int currentWindow = window;
        int scale = 2;
        if (var.ref.size()*scale > currentWindow) currentWindow = var.ref.size()*scale;
        for (const auto& a : var.alleles) {
            if (a.size() * scale > currentWindow) {
                currentWindow = a.size() * scale;
            }
        }

        // while the entropy of either flank is < some target entropy (~1 is fine), increase the flank sizes
        while (currentWindow < 2000) { // limit to one step > than this
            string refTarget = fastaReference.getSubSequence(var.sequenceName, var.position - 1 - currentWindow/2, currentWindow);
            if (entropy(refTarget.substr(0, refTarget.size()/2)) < 1 ||
                entropy(refTarget.substr(refTarget.size()/2)) < 1) {
                currentWindow *= scale;
            } else {
                break;
            }
        }

        // do the alignments
        getAlignment(var, fastaReference, ref, alignments, currentWindow);

        // stably left align the alignments
        for (vector<AltAlignment>::iterator a = alignments.begin(); a != alignments.end(); ++a) {
            Cigar cigarBefore = a->cigar;
            //cerr << a->seq << endl;
            //cerr << "before : " << a->pos << " " << joinCigar(a->cigar) << endl;
            long int prev = a->pos;
            stablyLeftAlign(a->seq, ref, a->cigar, 20, false);
            //cerr << "after  : " << a->pos << " " << joinCigar(a->cigar) << endl;
            if (a->pos != prev) cerr << "modified alignment @ " << var << endl;
        }
        //cout << var << endl;

        // transform the mappings
        // chop off leading matching bases
        // find the range of bp in the alleles
        // make the new ref allele
        // make the new alt alleles
        // emit the var

        long int newPosition = var.position+currentWindow/2;
        long int newEndPosition = var.position-currentWindow/2;
        // check for no-indel case
        int newLength = var.ref.size();
        bool giveUp = false;
        for (vector<AltAlignment>::iterator a = alignments.begin(); a != alignments.end() && !giveUp; ++a) {
            // get the first mismatching position
            Cigar::iterator c = a->cigar.begin();

            int rp = 0;
            int sp = 0;
            bool hitMismatch = false;

            int matchingBpAtStart = 0;
            int matchingBpAtEnd = 0;
            // will be set to true if the first reference position match is broken by a SNP, not an indel
            bool leadingSNP = false;

            while (c != a->cigar.end()) {
                char op = c->second;
                if (c == a->cigar.begin()) {
                    if (op != 'M') {
                        cerr << "alignment does not start on matched sequence" << endl;
                        cerr << var << endl;
                        exit(1);
                    }
                    int i = 0;
                    for ( ; i < c->first; ++i) {
                        if (ref[i] != a->seq[i]) {
                            leadingSNP = true;
                            break;
                        }
                    }
                    matchingBpAtStart = i;
                }
                if (!leadingSNP && c == (a->cigar.begin()+1)) {
                    // if the first thing we run into is an indel, step back, per VCF spec
                    if (op == 'D' || op == 'I') {
                        --matchingBpAtStart;
                    }
                }
                if (c == (a->cigar.end()-1)) {
                    if (op != 'M') {
                        // soft clip at end
                        // it'll be hard to interpret this
                        // the alignments sometimes generate this
                        // best thing to do is to move on
                        //cerr << "alignment does not end on matched sequence" << endl;
                        //cout << var << endl;
                        //exit(1);
                        giveUp = true;
                        break;
                    }
                    int i = 0;
                    for ( ; i < c->first; ++i) {
                        if (ref[ref.size()-1-i] != a->seq[a->seq.size()-1-i]) {
                            break;
                        }
                    }
                    matchingBpAtEnd = i;
                }
                ++c;
            }

            int altMismatchLength = a->seq.size() - matchingBpAtEnd - matchingBpAtStart;
            int refMismatchLength = (var.ref.size() + currentWindow) - matchingBpAtEnd - matchingBpAtStart;
            //cerr << "alt mismatch length " << altMismatchLength << endl
            //     << "ref mismatch length " << refMismatchLength << endl;
            long int newStart = var.position - currentWindow/2 + matchingBpAtStart;
            long int newEnd = newStart + refMismatchLength;
            //cerr << "ref should run from " << newStart << " to " << newStart + refMismatchLength << endl;
            newPosition = min(newStart, newPosition);
            newEndPosition = max(newEnd, newEndPosition);
            //cerr << newPosition << " " << newEndPosition << endl;
            //if (newRefSize < refMismatchLength) newRefSize = refMismatchLength;
        }

        // the alignment failed for some reason, continue
        if (giveUp) {
            cout << var << endl;
            continue;
        }

        //cerr << "new ref start " << newPosition << " and end " << newEndPosition << " was " << var.position << "," << var.position + var.ref.size() << endl;
        int newRefSize = newEndPosition - newPosition;
        string newRef = fastaReference.getSubSequence(var.sequenceName, newPosition-1, newRefSize);
        // get the number of bp to strip from the alts
        int stripFromStart = currentWindow/2 - (var.position - newPosition);
        int stripFromEnd = (currentWindow + newRefSize) - (stripFromStart + newRefSize) + (var.ref.size() - newRefSize);

        //cerr << "strip from start " << stripFromStart << endl;
        //cerr << "strip from end " << stripFromEnd << endl;

        vector<string> newAlt;
        vector<string>::iterator l = var.alt.begin();
        bool failedAlt = false;
        for (vector<AltAlignment>::iterator a = alignments.begin(); a != alignments.end();
             ++a, ++l) {
            int diff = newRef.size() - l->size();
            string alt = a->seq.substr(stripFromStart, a->seq.size() - (stripFromEnd + stripFromStart));
            newAlt.push_back(alt);
            if (alt.empty()) failedAlt = true;
        }

        // check the before/after haplotypes
        bool brokenRealignment = false;
        if (!newRef.empty() && !failedAlt) {
            int slop = 50; // 50 extra bp!
            int haplotypeStart = min(var.position, newPosition) - slop;
            int haplotypeEnd = max(var.position + var.ref.size(), newPosition + newRef.size()) + slop;
            string referenceHaplotype = fastaReference.getSubSequence(var.sequenceName, haplotypeStart - 1,
                                                                      haplotypeEnd - haplotypeStart);
            vector<string>::iterator o = var.alt.begin();
            vector<string>::iterator n = newAlt.begin();
            for ( ; o != var.alt.end() ; ++o, ++n) {
                // map the haplotypes
                string oldHaplotype = referenceHaplotype;
                string newHaplotype = referenceHaplotype;
                oldHaplotype.replace(var.position - haplotypeStart, var.ref.size(), *o);
                newHaplotype.replace(newPosition - haplotypeStart, newRef.size(), *n);
                if (oldHaplotype != newHaplotype) {
                    cerr << "broken left alignment!" << endl
                         << "old " << oldHaplotype << endl
                         << "new " << newHaplotype << endl;
                    cerr << "was: " << var << endl;
                    brokenRealignment = true;
                }
            }
        }

        // *if* everything is OK, update the variant
        if (!brokenRealignment && !newRef.empty() && !failedAlt) {
            var.ref = newRef;
            var.alt = newAlt;
            var.position = newPosition;
        }

        cout << var << endl;

        // for each parsedalternate, get the position
        // build a new vcf record for that position
        // unless we are already at the position !
        // take everything which is unique to that allele (records) and append it to the new record
        // then handle genotypes; determine the mapping between alleleic primitives and convert to phased haplotypes
        // this means taking all the parsedAlternates and, for each one, generating a pattern of allele indecies corresponding to it



        //for (vector<Variant>::iterator v = variants.begin(); v != variants.end(); ++v) {
    }

    return 0;

}
