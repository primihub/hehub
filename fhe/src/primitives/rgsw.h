/**
 * @file rgsw.h
 * @brief TODO
 */

#pragma once

#include "common/rnspolynomial.h"
#include "rlwe.h"
#include <array>

namespace hehub {

/**
 * @brief TODO
 *
 */
struct RgswSpec {};

using RgswCt = std::vector<RlweCt>;

/**
 * @brief TODO
 *
 * @param pt_ntt
 * @param sk
 * @param decomp_basis
 * @return RgswCt
 */
RgswCt rgsw_encrypt(const RlwePt &pt_ntt, const RlweSk &sk, 
                    const std::vector<std::vector<u64>> &decomp_basis);

/**
 * @brief TODO
 *
 * @param pt_ntt
 * @param sk
 * @param decomp_basis
 * @return RgswCt
 */
RgswCt rgsw_encrypt_montgomery(const RlwePt &pt_ntt, const RlweSk &sk, 
                               const std::vector<std::vector<u64>> &decomp_basis);

/**
 * @brief TODO
 *
 * @param pt
 * @param rgsw
 * @return RlweCt
 */
RlweCt ext_prod_montgomery(const RlwePt &pt, const RgswCt &rgsw);

} // namespace hehub
