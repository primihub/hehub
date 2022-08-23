/**
 * @file sampling.h
 * @brief Generate random polynomials.
 */

#include "rnspolynomial.h"

namespace hehub {

RnsPolynomial get_rand_ternary_poly(const size_t component_count,
                                    const size_t log_poly_len,
                                    const std::vector<u64> &moduli);

} // namespace hehub
