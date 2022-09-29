#include "ckks.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include <algorithm>
#include <numeric>

namespace hehub {

void rescale_by_one_prime_inplace(CkksCt &ct) {
    if (ct[0].modulus_vec() != ct[1].modulus_vec()) {
        throw std::invalid_argument(
            "Ill-formed ciphertext: modulus sets mismatch.");
    }
    if (ct[0].poly_len() != ct[1].poly_len()) {
        throw std::invalid_argument(
            "Ill-formed ciphertext: polynomial lengths mismatch.");
    }
    if (ct[0].component_count() != ct[1].component_count()) {
        throw std::invalid_argument(
            "Ill-formed ciphertext: component numbers mismatch.");
    }
    if (ct[0].component_count() == 1) {
        throw std::invalid_argument("Unable to drop the only one prime.");
    }

    const auto ct_moduli = ct[0].modulus_vec();
    const auto poly_len = ct[0].poly_len();
    const auto log_poly_len = ct[0].log_poly_len();
    const auto ct_mod_count = ct[0].component_count();
    const auto q_last = ct_moduli[ct_mod_count - 1];
    const auto half_q_last = q_last / 2;
    std::vector<u64> q_last_mod_qk, inv_q_last_mod_qk;
    for (size_t k = 0; k < ct_mod_count - 1; k++) {
        q_last_mod_qk.push_back(q_last % ct_moduli[k]); // need opt?
        inv_q_last_mod_qk.push_back(inverse_mod_prime(q_last, ct_moduli[k]));
    }

    for (auto &rns_poly : ct) {
        RnsPolynomial last_comp{poly_len, 1, std::vector{q_last}};
        last_comp[0] = rns_poly[ct_mod_count - 1];
        intt_negacyclic_inplace_lazy(last_comp);
        batched_strict_reduce(q_last, poly_len, last_comp[0].data());

        auto &last_comp_coeffs = last_comp[0];
        RnsPolynomial remainder_q_last(poly_len, ct_mod_count - 1, ct_moduli);
        for (size_t k = 0; k < remainder_q_last.component_count(); k++) {
            remainder_q_last[k] = last_comp_coeffs;

            /* Heuristically the ciphertext moduli are close in size, hence the
             * reduction step after copying can be skipped. */
            // batched_barrett_lazy(
            //     ct_moduli[k], poly_len, remainder_q_last[k].data());

            // The remainder needs to be of smallest possible abs value, i.e. in
            // [-q_last/2, q_last/2).
            for (size_t i = 0; i < poly_len; i++) {
                if (last_comp_coeffs[i] >= half_q_last) {
                    remainder_q_last[k][i] += ct_moduli[k] - q_last_mod_qk[k];
                }
            }
        }
        ntt_negacyclic_inplace_lazy(remainder_q_last);

        rns_poly.remove_components();
        rns_poly -= remainder_q_last;
        rns_poly *= inv_q_last_mod_qk;
    }
}

void ckks::rescale_inplace(CkksCt &ct, size_t dropping_primes) {
    if (dropping_primes == 1) {
        rescale_by_one_prime_inplace(ct);
    } else if (dropping_primes >= 2) {
        // TODO case
    } else {
        // TODO: throw error
    }
}

} // namespace hehub