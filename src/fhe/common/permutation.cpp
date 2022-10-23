#include "permutation.h"
#include "range/v3/view/zip.hpp"

using namespace std;
using namespace ranges::views;

namespace hehub {

const size_t GALOIS_GEN = 3;

const vector<u32> &root_index_factors() {
    struct RootIndexFactors : public vector<u32> {
        RootIndexFactors()
            : vector(1 << 17) // enables a param n <= 2^16
        {
            (*this)[0] = 1;
            for (size_t i = 1; i < this->size(); i++) {
                // modulo 2^32, consistent to modulo another smaller 2-power
                (*this)[i] = (*this)[i - 1] * GALOIS_GEN;
            }
        }
    };

    static RootIndexFactors global_root_index_factors;
    return global_root_index_factors;
}

RnsPolynomial cycle(const RnsPolynomial &poly_ntt, const size_t step) {
    if (poly_ntt.rep_form != PolyRepForm::value) {
        throw invalid_argument("poly_ntt is expected to be in NTT value form");
    }

    const auto len = poly_ntt.dimension();
    const auto loglen = poly_ntt.log_dimension();
    const auto components = poly_ntt.component_count();
    RnsPolynomial cycled(len, components, poly_ntt.modulus_vec());
    cycled.rep_form = PolyRepForm::value;

    auto mask = (1 << (loglen + 1)) - 1; // for fast modulo 2*len
    auto &root_indices = root_index_factors();
    auto index_factor = root_indices[step] & mask;
    for (size_t i = 0; i < len / 2; i++) {
        auto old_root_index = root_indices[i] & mask;
        auto from_position =
            __bit_rev_naive_16((old_root_index - 1) / 2, loglen);
        auto new_root_index = old_root_index * index_factor & mask;
        auto to_position = __bit_rev_naive_16((new_root_index - 1) / 2, loglen);

        for (size_t k = 0; k < components; k++) {
            cycled[k][to_position] = poly_ntt[k][from_position];
            cycled[k][len - 1 - to_position] =
                poly_ntt[k][len - 1 - from_position];
        }
    }

    return cycled;
}

RnsPolynomial involution(const RnsPolynomial &poly_ntt) {
    if (poly_ntt.rep_form != PolyRepForm::value) {
        throw invalid_argument("poly_ntt is expected to be in NTT value form");
    }

    const auto len = poly_ntt.dimension();
    const auto components = poly_ntt.component_count();
    RnsPolynomial involution(len, components, poly_ntt.modulus_vec());
    involution.rep_form = PolyRepForm::value;
    for (auto [new_component, old_component] : zip(involution, poly_ntt)) {
        for (size_t i = 0; i < len; i++) {
            new_component[i] = old_component[len - 1 - i];
        }
    }

    return involution;
}

} // namespace hehub
