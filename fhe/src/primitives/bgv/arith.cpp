#include "bgv.h"
#include "common/rns_transform.h"
#include "common/ntt.h"

namespace hehub {

RlweCt bgv::add(const RlweCt &ct1, const RlweCt &ct2) {
    return ::hehub::add(ct1, ct2);
}

RlweCt bgv::add_plain(const RlweCt &ct, const RlwePt &pt) {
    auto pt_under_ct_mod = rns_base_transform(pt, ct[0].modulus_vec());
    ntt_negacyclic_inplace_lazy(pt_under_ct_mod);
    return ::hehub::add_plain_core(ct, pt_under_ct_mod);
}

RlweCt bgv::sub(const RlweCt &ct1, const RlweCt &ct2) {
    return ::hehub::sub(ct1, ct2);
}

RlweCt bgv::sub_plain(const RlweCt &ct, const RlwePt &pt) {
    auto pt_under_ct_mod = rns_base_transform(pt, ct[0].modulus_vec());
    ntt_negacyclic_inplace_lazy(pt_under_ct_mod);
    return ::hehub::sub_plain_core(ct, pt_under_ct_mod);
}

RlweCt bgv::mult_plain(const RlweCt &ct, const RlwePt &pt) {
    auto pt_under_ct_mod = rns_base_transform(pt, ct[0].modulus_vec());
    ntt_negacyclic_inplace_lazy(pt_under_ct_mod);
    return ::hehub::mult_plain_core(ct, pt_under_ct_mod);
}

} // namespace hehub
