#include "rgsw.h"
#include "common/mod_arith.h"
#include "common/ntt.h"

using namespace std;

namespace hehub {

RlweCt ext_prod(const RlwePt &pt, const RgswCt &rgsw) {
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
    for (size_t i = 0; i < original_components; i++) {
        if (extended_moduli[i] != moduli[i]) {
            throw invalid_argument("Moduli mismatch.");
        }
    }

    const auto poly_len = pt.poly_len();
    const auto log_poly_len = pt.log_poly_len();
    for (auto &rlwe_sample : rgsw) {
        for (auto &poly : rlwe_sample) {
            if (poly.poly_len() != poly_len) {
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
    PolyDimensions extended_poly_dim{poly_len, extended_components,
                                     extended_moduli};
    for (auto &rns_poly : decomposed) {
        rns_poly = RnsPolynomial(extended_poly_dim);
    }

    // The components on the diagonal are reserved
    for (size_t k = 0; k < original_components; k++) {
        decomposed[k][k] = pt[k];
    }

    auto pt_intt(pt);
    intt_negacyclic_inplace_lazy(pt_intt);
    strict_reduce(pt_intt);

    // Copy the components not on the diagonal
    for (size_t rns_pl_idx = 0; rns_pl_idx < original_components;
         rns_pl_idx++) {
        for (int compo_idx = 0; compo_idx < extended_components; compo_idx++) {
            if (compo_idx == rns_pl_idx) {
                continue;
            }
            // Copy the "rns_pl_idx"-th component of pt_intt
            decomposed[rns_pl_idx][compo_idx] = pt_intt[rns_pl_idx];
            ntt_negacyclic_inplace_lazy(
                log_poly_len, extended_moduli[compo_idx],
                decomposed[rns_pl_idx][compo_idx].data());
        }
    }

    RlweCt ct_tilde{RnsPolynomial(extended_poly_dim),
                    RnsPolynomial(extended_poly_dim)};
    u128 temp_sum[poly_len];

    // Multiply the matrices with the RGSW
    for (auto half : {0, 1}) {
        // modulo original moduli part
        for (int k = 0; k < original_components; k++) {
            fill(temp_sum, temp_sum + poly_len, (u128)0);
            for (int rns_pl_idx = 0; rns_pl_idx < original_components;
                 rns_pl_idx++) {
                for (int i = 0; i < poly_len; i++) {
                    temp_sum[i] += (u128)decomposed[rns_pl_idx][k][i] *
                                   rgsw[rns_pl_idx][half][k][i];
                }
            }
            batched_montgomery_128_lazy(moduli[k], poly_len, temp_sum,
                                        ct_tilde[half][k].data());
        }

        // modulo new modulus part
        fill(temp_sum, temp_sum + poly_len, (u128)0);
        for (int rns_pl_idx = 0; rns_pl_idx < original_components;
             rns_pl_idx++) {
            for (int i = 0; i < poly_len; i++) {
                temp_sum[i] += (u128)(*decomposed[rns_pl_idx].last())[i] *
                               (*rgsw[rns_pl_idx][half].last())[i];
            }
        }
        batched_montgomery_128_lazy(*extended_moduli.crbegin(), poly_len,
                                    temp_sum, ct_tilde[half].last()->data());

        // Set as NTT value form
        ct_tilde[half].rep_form = PolyRepForm::value;
    }

    return ct_tilde;
}

} // namespace hehub
