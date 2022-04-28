#include "allele.hpp"

namespace vcflib {

using namespace std;

ostream& operator<<(ostream& out, VariantAllele& var) {
    out << var.position << " " << var.ref << " -> " << var.alt;
    return out;
}

VariantAllele operator+(const VariantAllele& a, const VariantAllele& b) {
    return VariantAllele(a.ref + b.ref, a.alt + b.alt, a.position);
}

bool operator<(const VariantAllele& a, const VariantAllele& b) {
    return a.repr < b.repr;
}

}
