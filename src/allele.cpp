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
    return std::tie(a.position, a.ref, a.alt) < std::tie(b.position, b.ref, b.alt);

}

bool operator==(const VariantAllele& a, const VariantAllele& b) {
    return a.ref == b.ref && a.alt == b.alt && a.position == b.position;
}

bool VariantAllele::is_pure_indel(void) {
    return ref.size() > 0 &&  alt == "" || alt.size() > 0 && ref == "";
}

void shift_mid_left(VariantAllele& a, VariantAllele& b) {
    if (!b.is_pure_indel()) {
        a.alt.append(b.alt.substr(0,1));
        a.ref.append(b.ref.substr(0,1));
        b.alt = b.alt.substr(1);
        b.ref = b.ref.substr(1);
        ++b.position;
    } else {

    }
}

void shift_mid_right(VariantAllele& a, VariantAllele& b) {
    if (!a.is_pure_indel()) {
        b.alt = a.alt.substr(a.alt.size()-1,1) + b.alt;
        b.ref = a.ref.substr(a.ref.size()-1,1) + b.ref;
        a.alt = a.alt.substr(0,a.alt.size()-1);
        a.ref = a.ref.substr(0,a.alt.size()-1);
        --b.position;
    } else {
    }
}

}
