#include "ckks.h"
#include "common/ntt.h"

namespace hehub {

const double EPS = std::pow(2.0, -50);

auto check_scaling_factor = [](const auto &in1, const auto &in2) {
    if (std::abs(in1.scaling_factor - in2.scaling_factor) > EPS) {
        throw std::invalid_argument("The scaling factors mismatch");
    }
};

CkksCt ckks::add(const CkksCt &ct1, const CkksCt &ct2) {
    check_scaling_factor(ct1, ct2);
    CkksCt sum_ct = ::hehub::add(ct1, ct2);
    sum_ct.scaling_factor = ct1.scaling_factor;
    return sum_ct;
}

CkksCt ckks::add_plain(const CkksCt &ct, const CkksPt &pt) {
    check_scaling_factor(ct, pt);
    auto pt_ntt(pt);
    ntt_negacyclic_inplace_lazy(pt_ntt);
    CkksCt sum_ct = ::hehub::add_plain_core(ct, pt_ntt);
    sum_ct.scaling_factor = ct.scaling_factor;
    return sum_ct;
}

CkksCt ckks::sub(const CkksCt &ct1, const CkksCt &ct2) {
    check_scaling_factor(ct1, ct2);
    CkksCt diff_ct = ::hehub::sub(ct1, ct2);
    diff_ct.scaling_factor = ct1.scaling_factor;
    return diff_ct;
}

CkksCt ckks::sub_plain(const CkksCt &ct, const CkksPt &pt) {
    check_scaling_factor(ct, pt);
    auto pt_ntt(pt);
    ntt_negacyclic_inplace_lazy(pt_ntt);
    CkksCt diff_ct = ::hehub::sub_plain_core(ct, pt_ntt);
    diff_ct.scaling_factor = ct.scaling_factor;
    return diff_ct;
}

CkksCt ckks::mult_plain(const CkksCt &ct, const CkksPt &pt) {
    auto pt_ntt(pt);
    ntt_negacyclic_inplace_lazy(pt_ntt);
    CkksCt prod_ct = ::hehub::mult_plain_core(ct, pt_ntt);
    prod_ct.scaling_factor = ct.scaling_factor * pt.scaling_factor;
    return prod_ct;
}

CkksQuadraticCt ckks::mult_low_level(const CkksCt &ct1, const CkksCt &ct2) {
    CkksQuadraticCt ct_prod;
    ct_prod[0] = ct1[0] * ct2[0];
    ct_prod[1] = ct1[0] * ct2[1] + ct1[1] * ct2[0];
    ct_prod[2] = ct1[1] * ct2[1];
    ct_prod.scaling_factor = ct1.scaling_factor * ct2.scaling_factor;
    return ct_prod;
}

CkksCt ckks::relinearize(const CkksQuadraticCt &ct, const RlweKsk &relin_key) {
    CkksCt ct_new = ext_prod_montgomery(ct[2], relin_key);
    rescale_inplace(ct_new); // this rescaling step shouldn't
                             // modify scaling factor
    ct_new.scaling_factor = ct.scaling_factor;

    ct_new[0] += ct[0];
    ct_new[1] += ct[1];
    return ct_new;
}

CkksCt ckks::conjugate(const CkksCt &ct, const RlweKsk &conj_key) {
    auto involution = RlweCt{involution(ct[0]), involution(ct[1])};
    CkksCt ct_conj = ext_prod_montgomery(involution[1], conj_key);
    ckks::rescale_inplace(ct_conj);
    ct_conj.scaling_factor = ct.scaling_factor; // the scaling factor
                                                // should remain
    ct_conj[0] += involution[0];
    return ct_conj;
}

CkksCt ckks::rotate(const CkksCt &ct, const RlweKsk &rot_key, const size_t step) {
    auto rotated = RlweCt{cycle(ct[0], step), cycle(ct[1], step)};
    CkksCt ct_rot = ext_prod_montgomery(rotated[1], rot_key);
    ckks::rescale_inplace(ct_rot);
    ct_rot.scaling_factor = ct.scaling_factor; // the scaling factor
                                               // should remain
    ct_rot[0] += rotated[0];
    return ct_rot;
}

} // namespace hehub
