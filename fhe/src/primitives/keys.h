/**
 * @file keys.h
 * @brief TODO
 *
 */

#pragma once

#include "rgsw.h"
#include "rlwe.h"

namespace hehub {

/**
 * @brief TODO
 *
 */
struct RlweKsk : public RgswCt {
    using RgswCt::RgswCt;

    /// @brief TODO
    /// @param rgsw
    RlweKsk(RgswCt &&rgsw) : RgswCt(std::move(rgsw)) {}

    /// @brief TODO
    /// @param sk_new
    /// @param sk_orig
    /// @param additional_mod
    RlweKsk(const RlweSk &sk_new, const RlweSk &sk_orig,
            const u64 additional_mod);
};

/**
 * @brief Generate and return a relinearization key. A relinearize key is an
 * RGSW encryption for sk^2 under sk, where sk is the secret key.
 * @param sk The secret key.
 * @param additional_mod The additional modulus for extending the RNS of the
 * relinearization key and the ciphertext to relinearize.
 * @return RlweKsk
 */
inline RlweKsk get_relin_key(const RlweSk &sk, const u64 additional_mod) {
    return RlweKsk(sk * sk, sk, additional_mod);
}

} // namespace hehub
