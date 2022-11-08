/**
 * @file bgv.h
 * @brief TODO
 *
 */

#pragma once
#include "fhe/primitives/keys.h"
#include "fhe/primitives/rlwe.h"

namespace hehub {
namespace bgv {

/**
 * @brief TODO
 *
 */
using BgvPt = RlwePt;

/**
 * @brief TODO
 *
 */
struct BgvCt : public RlweCt {
    using RlweCt::RlweCt;

    /// @brief TODO
    /// @param other
    BgvCt(RlweCt &&other) : RlweCt(std::move(other)) {}

    /// @brief TODO
    u64 plain_modulus = 1;
};

/**
 * @brief TODO
 * (TODO: merge this into BgvCt)
 */
struct BgvQuadraticCt : public std::array<RnsPolynomial, 3> {
    using std::array<RnsPolynomial, 3>::array;

    /// @brief TODO
    u64 plain_modulus = 1;
};


/**
 * @brief TODO
 *
 * @param modulus
 * @param data
 * @param slot_count
 * @return BgvPt
 */
BgvPt simd_encode(const std::vector<u64> &data, const u64 modulus,
                  size_t slot_count = 0);

/**
 * @brief TODO
 *
 * @param pt
 * @param data_size
 * @return std::vector<u64>
 */
std::vector<u64> simd_decode(const BgvPt &pt, size_t data_size = 0);

/**
 * @brief TODO
 *
 * @param sk
 * @param lifting_factor
 * @param components
 * @return RlweCt
 */
RlweCt get_rlwe_sample_lift_noise(const RlweSk &sk, const u64 lifting_factor,
                                  size_t components = 0);

/**
 * @brief TODO
 *
 * @param pt
 * @param rlwe_sk
 * @param ct_moduli
 * @return BgvCt
 */
BgvCt encrypt(const BgvPt &pt, const RlweSk &rlwe_sk,
              std::vector<u64> ct_moduli = std::vector<u64>{});

/**
 * @brief TODO
 *
 * @param ct
 * @param rlwe_sk
 * @return BgvPt
 */
BgvPt decrypt(const BgvCt &ct, const RlweSk &rlwe_sk);

/**
 * @brief TODO
 *
 * @param ct1
 * @param ct2
 * @return BgvCt
 */
BgvCt add(const BgvCt &ct1, const BgvCt &ct2);

/**
 * @brief TODO
 *
 * @param ct
 * @param pt
 * @return BgvCt
 */
BgvCt add_plain(const BgvCt &ct, const BgvPt &pt);

/**
 * @brief TODO
 *
 * @param ct1
 * @param ct2
 * @return BgvCt
 */
BgvCt sub(const BgvCt &ct1, const BgvCt &ct2);

/**
 * @brief TODO
 *
 * @param ct
 * @param pt
 * @return BgvCt
 */
BgvCt sub_plain(const BgvCt &ct, const BgvPt &pt);

/**
 * @brief TODO
 *
 * @param ct
 * @param pt
 * @return BgvCt
 */
BgvCt mult_plain(const BgvCt &ct, const BgvPt &pt);

/**
 * @brief TODO
 *
 * @param ct1
 * @param ct2
 * @return BgvQuadraticCt
 */
BgvQuadraticCt mult_low_level(const BgvCt &ct1, const BgvCt &ct2);

/**
 * @brief TODO
 *
 * @param ct
 * @param relin_key
 * @return BgvCt
 */
BgvCt relinearize(const BgvQuadraticCt &ct, const RlweKsk &relin_key);

/**
 * @brief TODO
 *
 * @param ct
 * @param dropping_primes
 */
void mod_switch_inplace(BgvCt &ct, size_t dropping_primes = 1);

} // namespace bgv
} // namespace hehub

// export types
namespace hehub
{
    using BgvPt = bgv::BgvPt;
    using BgvCt = bgv::BgvCt;
} // namespace hehub
