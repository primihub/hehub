#pragma once

#include "rnspolynomial.h"

namespace hehub {

/**
 * @brief TODO
 *
 * @param input_poly
 * @param new_moduli
 * @return RnsPolynomial
 */
RnsPolynomial rns_base_transform(RnsPolynomial input_poly,
                                 const std::vector<u64> &new_moduli);

} // namespace hehub
