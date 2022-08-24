#include "rlwe.h"
#include "common/sampling.h"

namespace hehub {

RlweSk::RlweSk(const PolyDimensions &poly_dim)
    : RnsPolynomial(get_rand_ternary_poly(poly_dim)) {}

} // namespace hehub
