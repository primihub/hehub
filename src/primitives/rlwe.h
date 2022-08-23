/**
 * @file rlwe.h
 * @brief Basics of the RLWE scheme.
 */

#include "common/rnspolynomial.h"
#include <array>

namespace hehub {

using RlwePt = RnsPolynomial;

using RlweCt = std::array<RnsPolynomial, 2>;

/**
 * @brief RLWE ternary secret key, which is a polynomial with coefficients
 * randomly sampled from {-1, 0, 1}. For convenience this is represented in RNS
 * and transformed to NTT form after sampling.
 */
class RlweSk : public RnsPolynomial {
public:
    RlweSk();

    /// Initialize the polynomial parameters and sample ternary coefficients.
    RlweSk(const size_t components, const size_t log_poly_len,
           const std::vector<u64> &moduli);
};

} // namespace hehub
