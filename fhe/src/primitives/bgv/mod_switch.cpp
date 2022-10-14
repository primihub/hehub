#include "bgv.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include <algorithm>
#include <numeric>

namespace hehub {

void mod_drop_one_prime_inplace(RlweCt &ct, u64 plain_modulus) {
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
    const auto inv_t_mod_q_last = inverse_mod_prime(plain_modulus, q_last);
    std::vector<u64> q_last_mod_qk, inv_q_last_mod_qk;
    for (size_t k = 0; k < ct_mod_count - 1; k++) {
        q_last_mod_qk.push_back(q_last % ct_moduli[k]); // need opt?
        inv_q_last_mod_qk.push_back(inverse_mod_prime(q_last, ct_moduli[k]));
    }

    for (auto &rns_poly : ct) {
        RnsPolynomial last_comp{poly_len, 1, std::vector{q_last}};
        last_comp[0] = rns_poly[ct_mod_count - 1];
        intt_negacyclic_inplace_lazy(last_comp);
        last_comp *= inv_t_mod_q_last;
        batched_strict_reduce(q_last, poly_len, last_comp[0].data());

        auto &coeffs_last_with_inv_t = last_comp[0];
        RnsPolynomial subtract_part(poly_len, ct_mod_count - 1, ct_moduli);
        for (size_t k = 0; k < subtract_part.component_count(); k++) {
            subtract_part[k] = coeffs_last_with_inv_t;

            /* Heuristically the ciphertext moduli are close in size, hence the
             * reduction step after copying can be skipped. */
            // batched_barrett_lazy(
            //     ct_moduli[k], poly_len, subtract_part[k].data());

            // The multiple of t needs to be of smallest possible abs value.
            for (size_t i = 0; i < poly_len; i++) {
                if (coeffs_last_with_inv_t[i] >= half_q_last) {
                    subtract_part[k][i] += ct_moduli[k] - q_last_mod_qk[k];
                }
            }
        }
        subtract_part *= plain_modulus;
        ntt_negacyclic_inplace_lazy(subtract_part);

        rns_poly.remove_components();
        rns_poly -= subtract_part;
        rns_poly *= inv_q_last_mod_qk;
        rns_poly *= (q_last % plain_modulus);
    }
}

void bgv::mod_switch_inplace(BgvCt &ct, size_t dropping_primes) {
    if (dropping_primes == 1) {
        mod_drop_one_prime_inplace(ct, ct.plain_modulus);
    } else if (dropping_primes >= 2) {
        // TODO case
        throw "under development";
    } else {
        throw std::invalid_argument(
            "The number of primes to be dropped is not positive.");
    }
}

} // namespace hehub
