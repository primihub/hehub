#include "bgv.h"
#include "common/ntt.h"
#include "common/rns_transform.h"

namespace hehub {

BgvCt bgv::add(const BgvCt &ct1, const BgvCt &ct2) {
    if (ct1.plain_modulus != ct2.plain_modulus) {
        throw std::invalid_argument("Plain moduli mismatch.");
    }
    BgvCt sum_ct = ::hehub::add(ct1, ct2);
    sum_ct.plain_modulus = ct1.plain_modulus;
    return sum_ct;
}

BgvCt bgv::add_plain(const BgvCt &ct, const BgvPt &pt) {
    if (pt.component_count() != 1 || pt.modulus_at(0) != ct.plain_modulus) {
        throw std::invalid_argument("plain moduli mismatch.");
    }
    auto pt_under_ct_mod = rns_base_transform(pt, ct[0].modulus_vec());
    ntt_negacyclic_inplace_lazy(pt_under_ct_mod);
    BgvCt sum_ct = ::hehub::add_plain_core(ct, pt_under_ct_mod);
    sum_ct.plain_modulus = ct.plain_modulus;
    return sum_ct;
}

BgvCt bgv::sub(const BgvCt &ct1, const BgvCt &ct2) {
    if (ct1.plain_modulus != ct2.plain_modulus) {
        throw std::invalid_argument("Plain moduli mismatch.");
    }
    BgvCt diff_ct = ::hehub::sub(ct1, ct2);
    diff_ct.plain_modulus = ct1.plain_modulus;
    return diff_ct;
}

BgvCt bgv::sub_plain(const BgvCt &ct, const BgvPt &pt) {
    if (pt.component_count() != 1 || pt.modulus_at(0) != ct.plain_modulus) {
        throw std::invalid_argument("plain moduli mismatch.");
    }
    auto pt_under_ct_mod = rns_base_transform(pt, ct[0].modulus_vec());
    ntt_negacyclic_inplace_lazy(pt_under_ct_mod);
    BgvCt diff_ct = ::hehub::sub_plain_core(ct, pt_under_ct_mod);
    diff_ct.plain_modulus = ct.plain_modulus;
    return diff_ct;
}

BgvCt bgv::mult_plain(const BgvCt &ct, const BgvPt &pt) {
    if (pt.component_count() != 1 || pt.modulus_at(0) != ct.plain_modulus) {
        throw std::invalid_argument("plain moduli mismatch.");
    }
    auto pt_under_ct_mod = rns_base_transform(pt, ct[0].modulus_vec());
    ntt_negacyclic_inplace_lazy(pt_under_ct_mod);
    BgvCt prod_ct = ::hehub::mult_plain_core(ct, pt_under_ct_mod);
    prod_ct.plain_modulus = ct.plain_modulus;
    return prod_ct;
}

BgvQuadraticCt bgv::mult_low_level(const BgvCt &ct1, const BgvCt &ct2) {
    if (ct1.plain_modulus != ct2.plain_modulus) {
        throw std::invalid_argument("Plain moduli mismatch.");
    }
    BgvQuadraticCt prod_ct;
    prod_ct[0] = ct1[0] * ct2[0];
    prod_ct[1] = ct1[0] * ct2[1] + ct1[1] * ct2[0];
    prod_ct[2] = ct1[1] * ct2[1];
    prod_ct.plain_modulus = ct1.plain_modulus;
    return prod_ct;
}

BgvCt bgv::relinearize(const BgvQuadraticCt &ct, const RlweKsk &relin_key) {
    BgvCt ct_new = ext_prod_montgomery(ct[2], relin_key);
    mod_switch_inplace(ct_new);

    ct_new[0] += ct[0];
    ct_new[1] += ct[1];
    ct_new.plain_modulus = ct.plain_modulus;
    return ct_new;
}

} // namespace hehub
