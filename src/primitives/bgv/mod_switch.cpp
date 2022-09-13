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
    const auto inv_t_mod_q_last_harvay = (u128(inv_t_mod_q_last) << 64) / q_last;
    std::vector<u64> q_last_mod_qk, inv_q_last_mod_qk, inv_q_last_mod_qk_harvey;
    for (size_t k = 0; k < ct_mod_count - 1; k++) {
        q_last_mod_qk.push_back(q_last % ct_moduli[k]); // can be opt?

        auto temp = inverse_mod_prime(q_last, ct_moduli[k]);
        inv_q_last_mod_qk.push_back(temp);
        inv_q_last_mod_qk_harvey.push_back(((u128)temp << 64) / ct_moduli[k]);
    }

    // In some practical cases plaintext modulus is greater than some or all
    // of ciphertext moduli. Hence plain_modulus needs to be reduced before
    // generation of harvay constants.
    std::vector<u64> plain_mod_reduced, plain_mod_harvey;
    for (const auto &ct_mod : ct_moduli) {
        auto curr_plain_mod_reduced = plain_modulus % ct_mod; // need opt(?)
        plain_mod_reduced.push_back(curr_plain_mod_reduced);
        plain_mod_harvey.push_back(((u128)curr_plain_mod_reduced << 64) /
                                   ct_mod);
    }

    for (auto &rns_poly : ct) {
        RnsPolynomial::ComponentData coeffs_last(rns_poly[ct_mod_count - 1]);
        intt_negacyclic_inplace_lazy(log_poly_len, q_last, coeffs_last.data());
        for (auto &coeff : coeffs_last) {
            coeff = mul_mod_harvey_lazy(q_last, coeff, inv_t_mod_q_last,
                                        inv_t_mod_q_last_harvay);
        }
        batched_strict_reduce(q_last, poly_len, coeffs_last.data());

        RnsPolynomial subtract_part(poly_len, ct_mod_count - 1, ct_moduli);
        for (size_t k = 0; k < subtract_part.component_count(); k++) {
            auto curr_mod = subtract_part.modulus_at(k);
            subtract_part[k] = coeffs_last;

            /* Heuristically the ciphertext moduli are close in size, hence the
             * reduction step after copying can be skipped. */
            // batched_barrett_lazy(curr_mod, poly_len, subtract_part[k].data());

            for (size_t i = 0; i < poly_len; i++) {
                if (coeffs_last[i] >= half_q_last) {
                    subtract_part[k][i] += ct_moduli[k] - q_last_mod_qk[k];
                }
                subtract_part[k][i] = mul_mod_harvey_lazy(
                    curr_mod, subtract_part[k][i], plain_mod_reduced[k],
                    plain_mod_harvey[k]);
            }
        }
        ntt_negacyclic_inplace_lazy(subtract_part);

        rns_poly.remove_components();
        rns_poly -= subtract_part;

        for (size_t k = 0; k < rns_poly.component_count(); k++) {
            for (auto &coeff : rns_poly[k]) {
                coeff = mul_mod_harvey_lazy(ct_moduli[k], coeff,
                                            inv_q_last_mod_qk[k],
                                            inv_q_last_mod_qk_harvey[k]);
            }
        }
    }
}

void bgv::mod_switch_inplace(RlweCt &ct, u64 plain_modulus,
                             size_t dropping_primes) {
    if (dropping_primes == 1) {
        mod_drop_one_prime_inplace(ct, plain_modulus);
    } else if (dropping_primes >= 2) {
        // TODO case
    } else {
        // TODO: throw error
    }
}

} // namespace hehub
