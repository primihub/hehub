#include "bgv.h"
#include "common/mod_arith.h"
#include "common/rns_transform.h"
#include <algorithm>
#include <cmath>

namespace hehub {

RlwePt bgv::simd_encode(const std::vector<u64> &data, const u64 modulus,
                        size_t slot_count) {
    for (auto datum : data) {
        if (datum >= modulus) {
            throw std::invalid_argument(
                "Data not being valid Z_p elements with p = " +
                std::to_string(modulus) + ".");
        }
    }

    if (slot_count == 0) {
        slot_count = 1 << (size_t)std::ceil(std::log2(data.size()));
    }
    auto data_size = data.size();
    if (data_size > slot_count) {
        throw std::invalid_argument("Cannot encode " +
                                    std::to_string(data_size) + " data into " +
                                    std::to_string(slot_count) + " slots.");
    }

    // Pack the data into plaintext slots.
    RlwePt pt(RlweParams{slot_count, 1, std::vector<u64>{modulus}});
    pt.rep_form = PolyRepForm::value;
    auto &pt_poly = pt[0];
    std::copy(data.begin(), data.end(), pt_poly.data());
    std::fill(pt_poly.begin() + data_size, pt_poly.end(), 0);

    // Transform the plaintext into coeff form which is required by encryption.
    intt_negacyclic_inplace_lazy(pt);

    return pt;
}

std::vector<u64> bgv::simd_decode(const RlwePt &pt, size_t data_size) {
    if (data_size == 0) {
        data_size = pt.dimension();
    }

    if (pt.component_count() != 1) {
        throw std::invalid_argument("Plaintext is in RNS with multiple moduli. "
                                    "Use big int version of decoding.");
    }

    auto pt_copy(pt);
    ntt_negacyclic_inplace_lazy(pt_copy);
    reduce_strict(pt_copy);
    std::vector<u64> data(pt_copy[0].begin(), pt_copy[0].end());
    data.resize(data_size);

    return data;
}

RlweCt bgv::get_rlwe_sample_lift_noise(const RlweSk &sk,
                                       const u64 lifting_factor,
                                       size_t components) {
    if (components == 0) {
        components = sk.component_count(); // actual argument
    }

    // Assume lifting_factor is already ensured to be coprime with RLWE modulus,
    // then a random ring element multiplied by it is still uniformly random.

    auto rlwe_sample = get_rlwe_sample(sk, components);
    for (auto &rns_poly : rlwe_sample) {
        rns_poly *= lifting_factor;
    }

    return rlwe_sample;
}

BgvCt bgv::encrypt(const RlwePt &pt, const RlweSk &rlwe_sk,
                   std::vector<u64> ct_moduli) {
    auto pt_modulus = pt.modulus_at(0);

    if (ct_moduli.empty()) {
        ct_moduli = rlwe_sk.modulus_vec();
    }
    if (std::find(ct_moduli.begin(), ct_moduli.end(), pt_modulus) !=
        ct_moduli.end()) {
        throw std::logic_error(
            "Plaintext modulus needs to be coprime with ciphertext modulus"
            " (i.e. cannot belong to ct_moduli).");
    }
    // Get a noise-lifted RLWE sample
    auto [c0, c1] =
        get_rlwe_sample_lift_noise(rlwe_sk, pt_modulus, ct_moduli.size());

    // Migrate the input plaintext to under needed ciphertext moduli
    auto pt_under_ct_mod = rns_base_transform(pt, ct_moduli);
    ntt_negacyclic_inplace_lazy(pt_under_ct_mod);

    c0 += pt_under_ct_mod;

    BgvCt ct = RlweCt{std::move(c0), std::move(c1)};
    ct.plain_modulus = pt_modulus;
    return ct;
}

BgvPt bgv::decrypt(const BgvCt &ct, const RlweSk &rlwe_sk) {
    // Apply RLWE decryption, obtaining the plaintext under ciphertext moduli
    // (and in coefficient form).
    auto pt_under_ct_mod = hehub::decrypt_core(ct, rlwe_sk);

    // Migrate the decrypted plaintext back to under original modulus
    auto pt =
        rns_base_transform(pt_under_ct_mod, std::vector{ct.plain_modulus});
    return pt;
}

} // namespace hehub
