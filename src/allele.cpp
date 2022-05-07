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

}
