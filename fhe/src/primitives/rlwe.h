/**
 * @file rlwe.h
 * @brief Basics of the RLWE scheme.
 */

#pragma once

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
struct RlweSk : public RnsPolynomial {
    using RnsPolynomial::RnsPolynomial;

    /// Create an empty (uninitialized) RLWE secret key.
    RlweSk();

    /// @brief TODO
    /// @param rns_poly
    RlweSk(RnsPolynomial &&rns_poly) : RnsPolynomial(rns_poly) {}

    /// @brief Initialize the polynomial parameters and sample ternary
    /// coefficients. For efficiency the secret key is stored in NTT form.
    /// @param poly_dim
    RlweSk(const PolyDimensions &poly_dim);
};

/**
 * @brief TODO
 * 
 * @param sk 
 * @param components 
 * @return RlweCt 
 */
RlweCt get_rlwe_sample(const RlweSk &sk, size_t components = 0);

/**
 * @brief TODO
 *
 * @param pt
 * @param sk
 * @return RlweCt
 */
RlweCt encrypt_core(const RlwePt &pt, const RlweSk &sk);

/**
 * @brief TODO
 *
 * @param ct
 * @param sk
 * @return RlwePt
 */
RlwePt decrypt_core(const RlweCt &ct, const RlweSk &sk);

/**
 * @brief TODO
 *
 * @param ct1
 * @param ct2
 * @return RlweCt
 */
RlweCt add(const RlweCt &ct1, const RlweCt &ct2);

/**
 * @brief TODO
 *
 * @param ct
 * @param pt
 * @return RlweCt
 */
RlweCt add_plain_core(const RlweCt &ct, const RlwePt &pt);

/**
 * @brief TODO
 *
 * @param ct1
 * @param ct2
 * @return RlweCt
 */
RlweCt sub(const RlweCt &ct1, const RlweCt &ct2);

/**
 * @brief TODO
 *
 * @param ct
 * @param pt
 * @return RlweCt
 */
RlweCt sub_plain_core(const RlweCt &ct, const RlwePt &pt);

/**
 * @brief TODO
 *
 * @param ct
 * @param pt
 * @return RlweCt
 */
RlweCt mult_plain_core(const RlweCt &ct, const RlwePt &pt);

} // namespace hehub
