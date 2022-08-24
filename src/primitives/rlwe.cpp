#include "rlwe.h"
#include "common/sampling.h"

namespace hehub {

RlweSk::RlweSk(const PolyDimensions &poly_dim)
    : RnsPolynomial(get_rand_ternary_poly(poly_dim)) {
    // For efficiency the secret key is stored in NTT form.
    ntt_negacyclic_inplace_lazy(*this);
}

} // namespace hehub
