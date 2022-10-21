#include "rgsw.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include "range/v3/view/zip.hpp"

using namespace std;
using namespace ranges::views;

namespace hehub {

RgswCt rgsw_encrypt(const RlwePt &pt_ntt, const RlweSk &sk,
                    const vector<vector<u64>> &decomp_basis) {
    if (pt_ntt.rep_form == PolyRepForm::coeff) {
        throw invalid_argument("Plaintext is expected in NTT form.");
    }

    auto sample_count = decomp_basis.size();
    RgswCt rgsw(sample_count);

    // generate the masks
    for (auto &rlwe_sample : rgsw) {
        rlwe_sample = get_rlwe_sample(sk);
    }

    // add pt multiplied by decomposition basis
    for (auto [rlwe_sample, basis_component] : zip(rgsw, decomp_basis)) {
        rlwe_sample[0] += pt_ntt * basis_component;
    }

    return rgsw;
}

RgswCt rgsw_encrypt_montgomery(const RlwePt &pt_ntt, const RlweSk &sk,
                               const vector<vector<u64>> &decomp_basis) {
    auto rgsw = rgsw_encrypt(pt_ntt, sk, decomp_basis);

    auto get_2to64_reduced = [](const u64 modulus) {
        const u64 _2to64_but_one = (u64)(-1LL) % modulus;
        return _2to64_but_one + 1;
    };

    auto moduli = rgsw[0][0].modulus_vec();
    std::vector<u64> mont_consts; // montgomery constants, which is 2^64 % q
    for (auto &modulus : moduli) {
        mont_consts.push_back(get_2to64_reduced(modulus));
    }

    for (auto &rlwe_sample : rgsw) {
        for (auto &poly : rlwe_sample) {
            poly *= mont_consts;
        }
    }

    return rgsw;
}

RlweCt ext_prod_montgomery(const RlwePt &pt, const RgswCt &rgsw) {
    const auto &moduli = pt.modulus_vec();
    if (rgsw.empty()) {
        throw invalid_argument("Empty RGSW ciphertext.");
    }

    auto extended_moduli = rgsw[0][0].modulus_vec();
    const auto original_components = pt.component_count();
    const auto extended_components = original_components + 1;
    if (extended_moduli.size() < extended_components) {
        throw invalid_argument("Invalid component number in RGSW ciphertext.");
    }
    extended_moduli.resize(extended_components);
    *extended_moduli.rbegin() = *rgsw[0][0].modulus_vec().crbegin();
    for (auto [modulus_rgsw, modulus_pt] : zip(extended_moduli, moduli)) {
        if (modulus_rgsw != modulus_pt) {
            throw invalid_argument("Moduli mismatch.");
        }
    }

    const auto dimension = pt.dimension();
    const auto log_dimension = pt.log_dimension();
    for (auto &rlwe_sample : rgsw) {
        for (auto &poly : rlwe_sample) {
            if (poly.dimension() != dimension) {
                throw invalid_argument("Polynomial lengths mismatch.");
            }
            if (poly.component_count() != extended_components ||
                poly.modulus_vec() != extended_moduli) {
                throw invalid_argument("Inconsistent RGSW ciphertext.");
            }
        }
    }

    // The decomposed pt, which forms the component matrix
    vector<RnsPolynomial> decomposed(original_components);
    RlweParams extended_params{dimension, extended_components, extended_moduli};
    for (auto &rns_poly : decomposed) {
        rns_poly = RnsPolynomial(extended_params);
    }

    // The components on the diagonal are reserved
    for (size_t k = 0; k < original_components; k++) {
        decomposed[k][k] = pt[k];
    }

    auto pt_intt(pt);
    intt_negacyclic_inplace_lazy(pt_intt);
    reduce_strict(pt_intt);

    // Copy the components not on the diagonal
    for (size_t poly_idx = 0; poly_idx < original_components; poly_idx++) {
        for (int compo_idx = 0; compo_idx < extended_components; compo_idx++) {
            if (compo_idx == poly_idx) {
                continue;
            }
            // Copy the "poly_idx"-th component of pt_intt
            decomposed[poly_idx][compo_idx] = pt_intt[poly_idx];
            ntt_negacyclic_inplace_lazy(log_dimension,
                                        extended_moduli[compo_idx],
                                        decomposed[poly_idx][compo_idx].data());
        }
    }

    RlweCt ct_tilde{RnsPolynomial(extended_params),
                    RnsPolynomial(extended_params)};
    u128 temp_sum[dimension];

    // Multiply the matrices with the RGSW
    for (auto half : {0, 1}) {
        // modulo original moduli part
        for (int k = 0; k < original_components; k++) {
            fill(temp_sum, temp_sum + dimension, (u128)0);
            for (int poly_idx = 0; poly_idx < original_components; poly_idx++) {
                for (int i = 0; i < dimension; i++) {
                    temp_sum[i] += (u128)decomposed[poly_idx][k][i] *
                                   rgsw[poly_idx][half][k][i];
                }
            }
            batched_montgomery_128_lazy(moduli[k], dimension, temp_sum,
                                        ct_tilde[half][k].data());
        }

        // modulo new modulus part
        fill(temp_sum, temp_sum + dimension, (u128)0);
        for (int poly_idx = 0; poly_idx < original_components; poly_idx++) {
            for (int i = 0; i < dimension; i++) {
                temp_sum[i] += (u128)(*decomposed[poly_idx].last())[i] *
                               (*rgsw[poly_idx][half].last())[i];
            }
        }
        batched_montgomery_128_lazy(*extended_moduli.crbegin(), dimension,
                                    temp_sum, ct_tilde[half].last()->data());

        // Set as NTT value form
        ct_tilde[half].rep_form = PolyRepForm::value;
    }

    return ct_tilde;
}

} // namespace hehub
