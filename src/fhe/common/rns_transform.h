#pragma once

#include "rns.h"

namespace hehub {

/**
 * @brief TODO
 * @note Currently this function accepts input of type RnsPolynomial, which can
 * be generalized to RnsIntVec if useful.
 * @param input_poly
 * @param new_moduli
 * @return RnsPolynomial
 */
RnsPolynomial rns_base_transform(RnsPolynomial input_poly,
                                 const std::vector<u64> &new_moduli);

} // namespace hehub
