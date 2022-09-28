#include "bigintpoly.h"

namespace hehub
{

UBigIntPoly::UBigIntPoly(const RnsPolynomial &rns_poly) {
    const auto poly_len(rns_poly.poly_len());
    const auto component_count(rns_poly.component_count());
    CRTComposer crt_composer(rns_poly.modulus_vec());
    for (size_t i = 0; i < poly_len; i++) {
        std::vector<u64> remainder_coeffs;
        for (size_t j = 0; j < component_count; j++) {
            remainder_coeffs.push_back(rns_poly[j][i]);
        }
        coeffs_.push_back(crt_composer.compose(remainder_coeffs));
    }
}

std::ostream &operator<<(std::ostream &out, const UBigIntPoly &big_int_poly) {
    for (size_t i = big_int_poly.coeffs_.size() - 1; i >= 0; i--) {
        out << big_int_poly.coeffs_[i];
        if (i == 0) { break; }
        out << "*X^" << i << " + ";
    }
    return out;
}

} // namespace hehub

