#include "rlwe.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include "common/sampling.h"

namespace hehub {

RlweSk::RlweSk(const PolyDimensions &poly_dim)
    : RnsPolynomial(get_rand_ternary_poly(poly_dim)) {
    // For efficiency the secret key is stored in NTT form.
    ntt_negacyclic_inplace_lazy(*this);
}

RlweCt get_rlwe_sample(const RlweSk &sk, const PolyDimensions &poly_dim) {
    auto ax = get_rand_uniform_poly(poly_dim, PolyRepForm::value);
    auto ex = get_rand_gaussian_poly(poly_dim);
    auto bx = ex - ax * sk;

    return RlweCt{std::move(ax), std::move(bx)};
}

RlweCt encrypt(const RlweSk &sk, const RlwePt &pt) {
    const auto poly_len = pt.poly_len();
    const auto &moduli = pt.moduli_vec();
    const auto components = moduli.size();
    PolyDimensions poly_dim{poly_len, components, moduli};

    if (pt.rep_form == PolyRepForm::value) {
        throw std::invalid_argument("Plaintext not in coeff representation.");
    }
    auto pt_ntt(pt);
    ntt_negacyclic_inplace_lazy(pt_ntt);

    auto [ax, bx] = get_rlwe_sample(sk, poly_dim);
    bx += pt_ntt;

    return RlweCt{std::move(ax), std::move(bx)};
}

RlwePt decrypt(const RlweSk &sk, const RlweCt &ct) {
    auto &[ax, bx] = ct;
    auto pt = bx + ax * sk;
    intt_negacyclic_inplace_lazy(pt);
    strict_reduce(pt);
    return pt;
}

} // namespace hehub
