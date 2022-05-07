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
    friend bool operator==(const VariantAllele& a, const VariantAllele& b);
public:
    string ref;
    string alt;
    long position;
    VariantAllele(string const & r, string const & a, long p)
        : ref(r), alt(a), position(p) { }
    bool is_pure_indel(void);
};

}
