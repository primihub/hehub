#include "permutation.h"

namespace hehub {

const size_t GALOIS_GEN = 3;

std::vector<u32> &root_index_factors() {
    struct RootIndexFactors : public std::vector<u32> {
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

const std::vector<size_t> &dlog_mod_2power_table(size_t loglen) {
    static std::vector<std::vector<size_t>> tables;
    
    if (tables.size() < loglen + 1) {
        tables.resize(loglen + 1);
    }
    if (!tables[loglen].empty()) {
        return tables[loglen];
    }

    auto &table = tables[loglen];
    auto len = 1 << loglen;
    auto doubled_len = 1 << (loglen + 1);
    table.resize(doubled_len);
    auto mask = doubled_len - 1; // for fast reduction
    size_t exp = 1;
    size_t dlog = 0;
    do {
        table[exp] = dlog;
        table[doubled_len - exp] = len - 1 - dlog;
        exp = exp * GALOIS_GEN & mask;
        dlog++;
    } while (exp != 1);

    return table;        
}

RnsPolynomial cycle(const RnsPolynomial &poly_ntt, const size_t step) {
    if (poly_ntt.rep_form != PolyRepForm::value) {
        throw std::invalid_argument(
            "poly_ntt is expected to be in NTT value form");
    }

    const auto len = poly_ntt.poly_len();
    const auto loglen = poly_ntt.log_poly_len();
    const auto components = poly_ntt.component_count();
    RnsPolynomial cycled(len, components, poly_ntt.modulus_vec());
    cycled.rep_form = PolyRepForm::value;

    auto mask = (1 << (loglen + 1)) - 1; // for fast modulo 2*len
    auto index_factor = root_index_factors()[step] & mask;
    for (size_t i = 0; i < len; i++) {
        auto orig_root_index = __bit_rev_naive_16(i, loglen) * 2 + 1;
        auto new_root_index = orig_root_index * index_factor & mask;
        auto new_index = __bit_rev_naive_16((new_root_index - 1) / 2, loglen);

        for (size_t k = 0; k < components; k++) {
            cycled[k][new_index] = poly_ntt[k][i];
        }
    }

    return cycled;
}

RnsPolynomial involute(const RnsPolynomial &poly_ntt) {
    if (poly_ntt.rep_form != PolyRepForm::value) {
        throw std::invalid_argument(
            "poly_ntt is expected to be in NTT value form");
    }

    const auto len = poly_ntt.poly_len();
    const auto components = poly_ntt.component_count();
    RnsPolynomial involuted(len, components, poly_ntt.modulus_vec());
    involuted.rep_form = PolyRepForm::value;
    for (size_t k = 0; k < components; k++) {
        for (size_t i = 0; i < len; i++) {
            involuted[k][i] = poly_ntt[k][len - 1 - i];
        }
    }

    return involuted;
}

} // namespace hehub
