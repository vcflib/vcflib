#include "Variant.h"
#include "convert.h"
#include "join.h"
#include "split.h"
#include "Fasta.h"
#include <set>
#include <vector>
#include <getopt.h>
#include <cmath>

using namespace std;
using namespace vcflib;


// Attempts to left-realign all the indels represented by the alignment cigar.
//
// This is done by shifting all indels as far left as they can go without
// mismatch, then merging neighboring indels of the same class.  leftAlign
// updates the alignment cigar with changes, and returns true if realignment
// changed the alignment cigar.
//
// To left-align, we move multi-base indels left by their own length as long as
// the preceding bases match the inserted or deleted sequence.  After this
// step, we handle multi-base homopolymer indels by shifting them one base to
// the left until they mismatch the reference.
//
// To merge neighboring indels, we iterate through the set of left-stabilized
// indels.  For each indel we add a new cigar element to the new cigar.  If a
// deletion follows a deletion, or an insertion occurs at the same place as
// another insertion, we merge the events by extending the previous cigar
// element.
//
// In practice, we must call this function until the alignment is stabilized.

#define VCFLEFTALIGN_DEBUG(msg) \
    if (false) { cerr << msg; }

class VCFIndelAllele {
    friend ostream& operator<<(ostream&, const VCFIndelAllele&);
    friend bool operator==(const VCFIndelAllele&, const VCFIndelAllele&);
    friend bool operator!=(const VCFIndelAllele&, const VCFIndelAllele&);
    friend bool operator<(const VCFIndelAllele&, const VCFIndelAllele&);
public:
    bool insertion;
    int length;
    int position;
    int readPosition;
    string sequence;

    bool homopolymer(void);

    VCFIndelAllele(bool i, int l, int p, int rp, string s)
        : insertion(i), length(l), position(p), readPosition(rp), sequence(s)
        { }
};

bool FBhomopolymer(string sequence);
ostream& operator<<(ostream& out, const VCFIndelAllele& indel);
bool operator==(const VCFIndelAllele& a, const VCFIndelAllele& b);
bool operator!=(const VCFIndelAllele& a, const VCFIndelAllele& b);
bool operator<(const VCFIndelAllele& a, const VCFIndelAllele& b);

bool VCFIndelAllele::homopolymer(void) {
    string::iterator s = sequence.begin();
    char c = *s++;
    while (s != sequence.end()) {
        if (c != *s++) return false;
    }
    return true;
}

bool FBhomopolymer(string sequence) {
    string::iterator s = sequence.begin();
    char c = *s++;
    while (s != sequence.end()) {
        if (c != *s++) return false;
    }
    return true;
}

ostream& operator<<(ostream& out, const VCFIndelAllele& indel) {
    string t = indel.insertion ? "i" : "d";
    out << t <<  ":" << indel.position << ":" << indel.readPosition << ":" << indel.sequence;
    return out;
}

bool operator==(const VCFIndelAllele& a, const VCFIndelAllele& b) {
    return (a.insertion == b.insertion
            && a.length == b.length
            && a.position == b.position
            && a.sequence == b.sequence);
}

bool operator!=(const VCFIndelAllele& a, const VCFIndelAllele& b) {
    return !(a==b);
}

bool operator<(const VCFIndelAllele& a, const VCFIndelAllele& b) {
    ostringstream as, bs;
    as << a;
    bs << b;
    return as.str() < bs.str();
}


class AltAlignment {
public:
    unsigned int pos;
    string seq;
    vector<pair<int, string> > cigar;
    AltAlignment(unsigned int& p,
                 string& s,
                 string& c) {
        pos = p;
        seq = s;
        cigar = splitCigar(c);
    }
};

double entropy(const string& st) {
    vector<char> stvec(st.begin(), st.end());
    set<char> alphabet(stvec.begin(), stvec.end());
    vector<double> freqs;
    for (set<char>::iterator c = alphabet.begin(); c != alphabet.end(); ++c) {
        int ctr = 0;
        for (vector<char>::iterator s = stvec.begin(); s != stvec.end(); ++s) {
            if (*s == *c) {
                ++ctr;
            }
        }
        freqs.push_back((double)ctr / (double)stvec.size());
    }
    double ent = 0;
    double ln2 = log(2);
    for (vector<double>::iterator f = freqs.begin(); f != freqs.end(); ++f) {
        ent += *f * log(*f)/ln2;
    }
    ent = -ent;
    return ent;
}

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
    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
        string alt = pad + leftFlank + *a + rightFlank + pad;
        CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
        unsigned int referencePos;
        string cigar;
        sw.Align(referencePos, cigar, ref, alt);
        alignments.push_back(AltAlignment(referencePos, alt, cigar));
    }
}


bool stablyLeftAlign(string& alternateSequence, string referenceSequence, int maxiterations = 50, bool debug = false);
int countMismatches(string& alternateSequence, string referenceSequence);

bool leftAlign(string& alternateSequence, Cigar& cigar, string& referenceSequence, bool debug = false) {

    int arsOffset = 0; // pointer to insertion point in aligned reference sequence
    string alignedReferenceSequence = referenceSequence;
    int aabOffset = 0;
    string alignmentAlignedBases = alternateSequence;

    // store information about the indels
    vector<VCFIndelAllele> indels;

    int rp = 0;  // read position, 0-based relative to read
    int sp = 0;  // sequence position

    string softBegin;
    string softEnd;

    stringstream cigar_before, cigar_after;
    for (vector<pair<int, string> >::const_iterator c = cigar.begin();
        c != cigar.end(); ++c) {
        unsigned int l = c->first;
        char t = c->second.at(0);

        cigar_before << l << t;
        if (t == 'M') { // match or mismatch
            sp += l;
            rp += l;
        } else if (t == 'D') { // deletion
            indels.push_back(VCFIndelAllele(false, l, sp, rp, referenceSequence.substr(sp, l)));
            alignmentAlignedBases.insert(rp + aabOffset, string(l, '-'));
            aabOffset += l;
            sp += l;  // update reference sequence position
        } else if (t == 'I') { // insertion
            indels.push_back(VCFIndelAllele(true, l, sp, rp, alternateSequence.substr(rp, l)));
            alignedReferenceSequence.insert(sp + softBegin.size() + arsOffset, string(l, '-'));
            arsOffset += l;
            rp += l;
        } else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
            // remove these bases from the refseq and read seq, but don't modify the alignment sequence
            if (rp == 0) {
                alignedReferenceSequence = string(l, '*') + alignedReferenceSequence;
                softBegin = alignmentAlignedBases.substr(0, l);
            } else {
                alignedReferenceSequence = alignedReferenceSequence + string(l, '*');
                softEnd = alignmentAlignedBases.substr(alignmentAlignedBases.size() - l, l);
            }
            rp += l;
        } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
        } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
            sp += l;
        }
    }


    int alignedLength = sp;

    VCFLEFTALIGN_DEBUG("| " << cigar_before.str() << endl
       << "| " << alignedReferenceSequence << endl
       << "| " << alignmentAlignedBases << endl);

    // if no indels, return the alignment
    if (indels.empty()) { return false; }

    // for each indel, from left to right
    //     while the indel sequence repeated to the left and we're not matched up with the left-previous indel
    //         move the indel left

    vector<VCFIndelAllele>::iterator previous = indels.begin();
    for (vector<VCFIndelAllele>::iterator id = indels.begin(); id != indels.end(); ++id) {

        // left shift by repeats
        //
        // from 1 base to the length of the indel, attempt to shift left
        // if the move would cause no change in alignment optimality (no
        // introduction of mismatches, and by definition no change in gap
        // length), move to the new position.
        // in practice this moves the indel left when we reach the size of
        // the repeat unit.
        //
        int steppos, readsteppos;
        VCFIndelAllele& indel = *id;
        int i = 1;
        while (i <= indel.length) {

            int steppos = indel.position - i;
            int readsteppos = indel.readPosition - i;

#ifdef VERBOSE_DEBUG
            if (debug) {
                if (steppos >= 0 && readsteppos >= 0) {
                    cerr << referenceSequence.substr(steppos, indel.length) << endl;
                    cerr << alternateSequence.substr(readsteppos, indel.length) << endl;
                    cerr << indel.sequence << endl;
                }
            }
#endif
            while (steppos >= 0 && readsteppos >= 0
                   && indel.sequence == referenceSequence.substr(steppos, indel.length)
                   && indel.sequence == alternateSequence.substr(readsteppos, indel.length)
                   && (id == indels.begin()
                       || (previous->insertion && steppos >= previous->position)
                       || (!previous->insertion && steppos >= previous->position + previous->length))) {
                VCFLEFTALIGN_DEBUG((indel.insertion ? "insertion " : "deletion ") << indel << " shifting " << i << "bp left" << endl);
                indel.position -= i;
                indel.readPosition -= i;
                steppos = indel.position - i;
                readsteppos = indel.readPosition - i;
            }
            do {
                ++i;
            } while (i <= indel.length && indel.length % i != 0);
        }

        // left shift indels with exchangeable flanking sequence
        //
        // for example:
        //
        //    GTTACGTT           GTTACGTT
        //    GT-----T   ---->   G-----TT
        //
        // GTGTGACGTGT           GTGTGACGTGT
        // GTGTG-----T   ---->   GTG-----TGT
        //
        // GTGTG-----T           GTG-----TGT
        // GTGTGACGTGT   ---->   GTGTGACGTGT
        //
        //
        steppos = indel.position - 1;
        readsteppos = indel.readPosition - 1;
        while (steppos >= 0 && readsteppos >= 0
               && alternateSequence.at(readsteppos) == referenceSequence.at(steppos)
               && alternateSequence.at(readsteppos) == indel.sequence.at(indel.sequence.size() - 1)
               && (id == indels.begin()
                   || (previous->insertion && indel.position - 1 >= previous->position)
                   || (!previous->insertion && indel.position - 1 >= previous->position + previous->length))) {
            VCFLEFTALIGN_DEBUG((indel.insertion ? "insertion " : "deletion ") << indel << " exchanging bases " << 1 << "bp left" << endl);
            indel.sequence = indel.sequence.at(indel.sequence.size() - 1) + indel.sequence.substr(0, indel.sequence.size() - 1);
            indel.position -= 1;
            indel.readPosition -= 1;
            steppos = indel.position - 1;
            readsteppos = indel.readPosition - 1;
        }
        // tracks previous indel, so we don't run into it with the next shift
        previous = id;
    }

    // bring together floating indels
    // from left to right
    // check if we could merge with the next indel
    // if so, adjust so that we will merge in the next step
    if (indels.size() > 1) {
        previous = indels.begin();
        for (vector<VCFIndelAllele>::iterator id = (indels.begin() + 1); id != indels.end(); ++id) {
            VCFIndelAllele& indel = *id;
            // parsimony: could we shift right and merge with the previous indel?
            // if so, do it
            int prev_end_ref = previous->insertion ? previous->position : previous->position + previous->length;
            int prev_end_read = !previous->insertion ? previous->readPosition : previous->readPosition + previous->length;
            if (previous->insertion == indel.insertion
                    && ((previous->insertion
                        && (previous->position < indel.position
                        && previous->readPosition + previous->readPosition < indel.readPosition))
                        ||
                        (!previous->insertion
                        && (previous->position + previous->length < indel.position)
                        && (previous->readPosition < indel.readPosition)
                        ))) {
                if (previous->homopolymer()) {
                    string seq = referenceSequence.substr(prev_end_ref, indel.position - prev_end_ref);
                    string readseq = alternateSequence.substr(prev_end_read, indel.position - prev_end_ref);
                    VCFLEFTALIGN_DEBUG("seq: " << seq << endl << "readseq: " << readseq << endl);
                    if (previous->sequence.at(0) == seq.at(0)
                            && FBhomopolymer(seq)
                            && FBhomopolymer(readseq)) {
                        VCFLEFTALIGN_DEBUG("moving " << *previous << " right to " 
                                << (indel.insertion ? indel.position : indel.position - previous->length) << endl);
                        previous->position = indel.insertion ? indel.position : indel.position - previous->length;
                    }
                } 
                else {
                    int pos = previous->position;
                    while (pos < (int) referenceSequence.length() &&
                            ((previous->insertion && pos + previous->length <= indel.position)
                            ||
                            (!previous->insertion && pos + previous->length < indel.position))
                            && previous->sequence 
                                == referenceSequence.substr(pos + previous->length, previous->length)) {
                        pos += previous->length;
                    }
                    if (pos < previous->position &&
                        ((previous->insertion && pos + previous->length == indel.position)
                        ||
                        (!previous->insertion && pos == indel.position - previous->length))
                       ) {
                        VCFLEFTALIGN_DEBUG("right-merging tandem repeat: moving " << *previous << " right to " << pos << endl);
                        previous->position = pos;
                    }
                }
            }
            previous = id;
        }
    }

    // for each indel
    //     if ( we're matched up to the previous insertion (or deletion) 
    //          and it's also an insertion or deletion )
    //         merge the indels
    //
    // and simultaneously reconstruct the cigar

    Cigar newCigar;

    if (!softBegin.empty()) {
        newCigar.push_back(make_pair(softBegin.size(), "S"));
    }

    vector<VCFIndelAllele>::iterator id = indels.begin();
    VCFIndelAllele last = *id++;
    if (last.position > 0) {
        newCigar.push_back(make_pair(last.position, "M"));
        newCigar.push_back(make_pair(last.length, (last.insertion ? "I" : "D")));
    } else {
        newCigar.push_back(make_pair(last.length, (last.insertion ? "I" : "D")));
    }
    int lastend = last.insertion ? last.position : (last.position + last.length);
    VCFLEFTALIGN_DEBUG(last << ",");

    for (; id != indels.end(); ++id) {
        VCFIndelAllele& indel = *id;
        VCFLEFTALIGN_DEBUG(indel << ",");
        if (indel.position < lastend) {
            cerr << "impossibility?: indel realigned left of another indel" << endl
                 << referenceSequence << endl << alternateSequence << endl;
            exit(1);
        } else if (indel.position == lastend && indel.insertion == last.insertion) {
            pair<int, string>& op = newCigar.back();
            op.first += indel.length;
        } else if (indel.position >= lastend) {  // also catches differential indels, but with the same position
            newCigar.push_back(make_pair(indel.position - lastend, "M"));
            newCigar.push_back(make_pair(indel.length, (indel.insertion ? "I" : "D")));
        }
        last = *id;
        lastend = last.insertion ? last.position : (last.position + last.length);
    }
    
    if (lastend < alignedLength) {
        newCigar.push_back(make_pair(alignedLength - lastend, "M"));
    }

    if (!softEnd.empty()) {
        newCigar.push_back(make_pair(softEnd.size(), "S"));
    }

    VCFLEFTALIGN_DEBUG(endl);

    cigar = newCigar;

    for (vector<pair<int, string> >::const_iterator c = cigar.begin();
        c != cigar.end(); ++c) {
        unsigned int l = c->first;
        char t = c->second.at(0);
        cigar_after << l << t;
    }

    //cerr << cigar_before.str() << " changes to " << cigar_after.str() << endl;
    VCFLEFTALIGN_DEBUG(cigar_after.str() << endl);

    // check if we're realigned
    if (cigar_after.str() == cigar_before.str()) {
        return false;
    } else {
        return true;
    }

}

// Iteratively left-aligns the indels in the alignment until we have a stable
// realignment.  Returns true on realignment success or non-realignment.
// Returns false if we exceed the maximum number of realignment iterations.
//
bool stablyLeftAlign(string& alternateSequence, string referenceSequence, Cigar& cigar, int maxiterations, bool debug) {

    if (!leftAlign(alternateSequence, cigar, referenceSequence, debug)) {

        return true;

    } else {

        bool result = true;
        while ((result = leftAlign(alternateSequence, cigar, referenceSequence, debug)) && --maxiterations > 0) { 
        }

        if (maxiterations <= 0) {
            return false;
        } else {
            return true;
        }

    }

}


void printSummary(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [file]" << endl
         << endl
         << "options:" << endl
         << "    -r, --reference FILE  Use this reference as a basis for realignment." << endl
         << "    -w, --window N        Use a window of this many bp when left aligning (150)." << endl
         << endl
         << "Left-aligns variants in the specified input file or stdin.  Window size is determined" << endl
         << "dynamically according to the entropy of the regions flanking the indel.  These must have" << endl
         << "entropy > 1 bit/bp, or be shorter than ~5kb." << endl;
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
        for (vector<string>::iterator a = var.alleles.begin(); a != var.alleles.end(); ++a) {
            if (a->size()*scale > currentWindow) {
                currentWindow = a->size()*scale;
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
                char op = c->second[0];
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

