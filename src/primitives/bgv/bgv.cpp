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
    RlwePt pt(PolyDimensions{slot_count, 1, std::vector<u64>{modulus}});
    pt.rep_form = PolyRepForm::value;
    auto mod_ptr = pt.moduli_vec().begin();
    auto &pt_poly = pt[0];
    std::copy(data.begin(), data.end(), pt_poly.data());
    std::fill(pt_poly.begin() + data_size, pt_poly.end(), 0);

    // Transform the plaintext into coeff form which is required by encryption.
    intt_negacyclic_inplace_lazy(pt);

    return pt;
}

std::vector<u64> bgv::simd_decode(RlwePt pt, size_t data_size) {
    if (data_size == 0) {
        data_size = pt.poly_len();
    }

    if (pt.component_count() != 1) {
        throw std::invalid_argument("Plaintext is in RNS with multiple moduli. "
                                    "Use big int version of decoding.");
    }

    ntt_negacyclic_inplace_lazy(pt);
    strict_reduce(pt);
    std::vector<u64> data(pt[0].begin(), pt[0].end());
    data.resize(data_size);

    return data;
}

RlweCt bgv::get_rlwe_sample_lift_noise(const RlweSk &sk,
                                       const PolyDimensions &poly_dim,
                                       const u64 lifting_factor) {
    auto rlwe_sample = get_rlwe_sample(sk, poly_dim);
    for (auto &rns_poly : rlwe_sample) {
        auto mod_ptr = rns_poly.moduli_vec().begin();
        for (auto &component : rns_poly) {
            auto component_mod = *(mod_ptr++);
            u64 lifting_factor_reduced = lifting_factor % component_mod;
            u64 lifting_factor_harvey =
                ((u128)lifting_factor_reduced << 64) / component_mod;
            for (auto &value : component) {
                value = mul_mod_harvey_lazy(component_mod, value,
                                            lifting_factor_reduced,
                                            lifting_factor_harvey);
            }
        }
    }

    return rlwe_sample;
}

RlweCt bgv::encrypt(const RlwePt &pt, const RlweSk &rlwe_sk,
                    std::vector<u64> ct_moduli) {
    if (pt.poly_len() != rlwe_sk.poly_len()) {
        throw std::invalid_argument(
            "The length of BGV plaintext and secret key mismatch. \nNote: "
            "may try to specify slot_count when encoding.");
    }
    auto pt_modulus = pt.modulus_at(0);

    if (ct_moduli.empty()) { ct_moduli = rlwe_sk.moduli_vec(); }
    if (std::find(ct_moduli.begin(), ct_moduli.end(), pt_modulus) !=
        ct_moduli.end()) {
        throw std::logic_error(
            "Plaintext modulus needs to be coprime with ciphertext modulus"
            " (i.e. cannot belong to ct_moduli).");
    }
    PolyDimensions poly_dim{pt.poly_len(), ct_moduli.size(), ct_moduli};
    // Get a noise-lifted RLWE sample
    auto [ax, bx] = get_rlwe_sample_lift_noise(rlwe_sk, poly_dim, pt_modulus);

    // Migrate the input plaintext to under needed ciphertext moduli
    auto pt_under_ct_mod = rns_base_transform(pt, ct_moduli);
    ntt_negacyclic_inplace_lazy(pt_under_ct_mod);

    bx += pt_under_ct_mod;
    return RlweCt{std::move(ax), std::move(bx)};
}

RlwePt bgv::decrypt(const RlweCt &ct, const RlweSk &rlwe_sk, u64 pt_modulus) {
    // Apply RLWE decryption, obtaining the plaintext under ciphertext moduli
    // (and in coefficient form).
    auto pt_under_ct_mod = hehub::decrypt(ct, rlwe_sk);

    // Migrate the decrypted plaintext back to under original modulus
    auto pt = rns_base_transform(pt_under_ct_mod, std::vector{pt_modulus});
    return pt;
}

} // namespace hehub
