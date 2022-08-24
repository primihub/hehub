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
    /// For efficiency the secret key is stored in NTT form.
    RlweSk(const PolyDimensions &poly_dim);
};

/**
 * @brief TODO
 * 
 * @param sk 
 * @param poly_dim 
 * @return RlweCt 
 */
RlweCt get_rlwe_sample(const RlweSk &sk, const PolyDimensions &poly_dim);

/**
 * @brief TODO
 *
 * @param sk
 * @param pt
 * @return RlweCt
 */
RlweCt encrypt(const RlweSk &sk, const RlwePt &pt);

/**
 * @brief TODO
 *
 * @param sk
 * @param ct
 * @return RlwePt
 */
RlwePt decrypt(const RlweSk &sk, const RlweCt &ct);

} // namespace hehub
