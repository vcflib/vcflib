#pragma once

#include <iostream>
#include <string>
#include <sstream>

namespace vcflib {

using namespace std;

/*
    A VariantAllele simply tracks [position,ref,alt] and has a string representation in 'repr'
*/

class VariantAllele {
    friend ostream& operator<<(ostream& out, VariantAllele& var);
    friend bool operator<(const VariantAllele& a, const VariantAllele& b);
    friend VariantAllele operator+(const VariantAllele& a, const VariantAllele& b);
public:
    string ref;
    string alt;
    string repr;
    long position;
    /* // TODO
    bool isSNP(void);
    bool isMNP(void);
    bool isInsertion(void);
    bool isDeletion(void);
    bool isIndel(void);
    */
    VariantAllele(string const & r, string const & a, long p)
        : ref(r), alt(a), position(p)
    {
        stringstream s;
        s << position << ":" << ref << "/" << alt;
        repr = s.str();
    }
};

}
