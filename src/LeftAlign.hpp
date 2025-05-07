#pragma once

#include "cigar.hpp"

#include <ostream>
#include <vector>

namespace vcflib {

using namespace std;

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
        : insertion(i), length(l), position(p), readPosition(rp), sequence(std::move(s))
        { }
};

bool FBhomopolymer(string sequence);
ostream& operator<<(ostream& out, const VCFIndelAllele& indel);
bool operator==(const VCFIndelAllele& a, const VCFIndelAllele& b);
bool operator!=(const VCFIndelAllele& a, const VCFIndelAllele& b);
bool operator<(const VCFIndelAllele& a, const VCFIndelAllele& b);

class AltAlignment {
public:
    unsigned int pos;
    string seq;
    vector<pair<int, char> > cigar;
    AltAlignment(unsigned int& p,
                 string& s,
                 string& c) {
        pos = p;
        seq = s;
        cigar = splitCigar(c);
    }
};

double entropy(const string& st);

int countMismatches(string& alternateSequence, string referenceSequence);

bool leftAlign(const string& alternateSequence, Cigar& cigar, string& referenceSequence, bool debug);

bool stablyLeftAlign(string& alternateSequence, string referenceSequence, Cigar& cigar, int maxiterations, bool debug);


}
