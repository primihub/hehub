#include "rlwe.h"
#include "common/sampling.h"

namespace hehub {

RlweSk::RlweSk(const size_t components, const size_t log_poly_len,
               const std::vector<u64> &moduli)
    : RnsPolynomial(get_rand_ternary_poly(components, log_poly_len, moduli)) {}

} // namespace hehub
