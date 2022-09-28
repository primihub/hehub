/**
 * @file rgsw.h
 * @brief TODO
 */

#include "common/rnspolynomial.h"
#include "rlwe.h"
#include <array>

namespace hehub {

/**
 * @brief TODO
 * 
 */
struct RgswSpec {

};

using RgswCt = std::vector<RlweCt>;

/**
 * @brief TODO
 * 
 * @param sk 
 * @param poly_dim 
 * @return RlweCt 
 */
RlweCt get_rgsw_sample(const RlweSk &sk, const PolyDimensions &poly_dim);

/**
 * @brief TODO
 *
 * @param sk
 * @param pt
 * @return RlweCt
 */
RlweCt rgsw_encrypt(const RlweSk &sk, const RlwePt &pt);

} // namespace hehub
