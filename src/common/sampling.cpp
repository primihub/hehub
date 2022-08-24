#include "sampling.h"
#include "rnspolynomial.h"
#include "type_defs.h"
#include <random>

namespace hehub {

std::random_device rand_dvc;
std::default_random_engine rand_engine;
std::uniform_int_distribution rand_ternary((i8)-1, (i8)1);

RnsPolynomial get_rand_ternary_poly(const PolyDimensions &poly_dim) {
    RnsPolynomial poly(poly_dim);
    auto poly_len = poly_dim.poly_len;
    auto component_count = poly_dim.component_count;

    // Sampling.
    i8 temp[poly_len];
    for (size_t i = 0; i < poly_len; i++) {
        temp[i] = rand_ternary(rand_engine);
    }

    // Transform to RNS representation.
    auto moduli_minus_one(poly.moduli_vec());
    for (auto &elem : moduli_minus_one) {
        elem--;
    }
    for (size_t i = 0; i < poly_len; i++) {
        if (temp[i] == -1) {
            for (size_t j = 0; j < component_count; j++) {
                poly[j][i] = moduli_minus_one[j];
            }
        } else {
            for (auto &component : poly.components()) {
                component[i] = temp[i];
            }
        }
    }

    return poly;
}

} // namespace hehub
