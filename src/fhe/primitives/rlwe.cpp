#include "rlwe.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include "common/sampling.h"

namespace hehub {

RlweSk::RlweSk(const RlweParams &params)
    : RnsPolynomial(get_rand_ternary_poly(params)) {}

RlweCt get_rlwe_sample(const RlweSk &sk, size_t components) {
    if (components == 0) {
        components = sk.component_count(); // the actual argument
    }
    RlweParams params{sk.dimension(), components, sk.modulus_vec()};
#ifdef HEHUB_DEBUG_RLWE_ZERO_C1
    auto c1 = get_zero_poly(params);
#else
    auto c1 = get_rand_uniform_poly(params, PolyRepForm::value);
#endif

#ifdef HEHUB_DEBUG_RLWE_ZERO_E
    auto ex = get_zero_poly(params);
#else
    auto ex = get_rand_gaussian_poly(params);
#endif

    auto c0 = ex - c1 * sk;
    return RlweCt{std::move(c0), std::move(c1)};
}

RlweCt encrypt_core(const RlwePt &pt, const RlweSk &sk) {
    const auto dimension = pt.dimension();
    const auto &moduli = pt.modulus_vec();
    const auto components = moduli.size();
    RlweParams params{dimension, components, moduli};

    if (pt.rep_form == PolyRepForm::value) {
        throw std::invalid_argument("Plaintext not in coeff representation.");
    }
    auto pt_ntt(pt);
    ntt_negacyclic_inplace_lazy(pt_ntt);

    auto [c0, c1] = get_rlwe_sample(sk, components);
    c0 += pt_ntt;

    return RlweCt{std::move(c0), std::move(c1)};
}

RlwePt decrypt_core(const RlweCt &ct, const RlweSk &sk) {
    auto &[c0, c1] = ct;
    auto pt = c0 + c1 * sk;

    // the obtained plaintext is now in NTT value representation
    intt_negacyclic_inplace_lazy(pt);
    strict_reduce(pt);
    return pt;
}

RlweCt add(const RlweCt &ct1, const RlweCt &ct2) {
    return RlweCt{ct1[0] + ct2[0], ct1[1] + ct2[1]};
}

RlweCt add_plain_core(const RlweCt &ct, const RlwePt &pt) {
    return RlweCt{ct[0] + pt, ct[1]};
}

RlweCt sub(const RlweCt &ct1, const RlweCt &ct2) {
    return RlweCt{ct1[0] - ct2[0], ct1[1] - ct2[1]};
}

RlweCt sub_plain_core(const RlweCt &ct, const RlwePt &pt) {
    return RlweCt{ct[0] - pt, ct[1]};
}

RlweCt mult_plain_core(const RlweCt &ct, const RlwePt &pt) {
    return RlweCt{ct[0] * pt, ct[1] * pt};
}

} // namespace hehub
