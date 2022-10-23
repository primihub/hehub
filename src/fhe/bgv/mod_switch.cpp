#include "bgv.h"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include "range/v3/view/zip.hpp"
#include <algorithm>
#include <numeric>

using namespace ranges::views;

namespace hehub {

void mod_drop_one_prime_inplace(RlweCt &ct, u64 plain_modulus) {
    if (ct[0].modulus_vec() != ct[1].modulus_vec()) {
        throw std::invalid_argument(
            "Ill-formed ciphertext: modulus sets mismatch.");
    }
    if (ct[0].dimension() != ct[1].dimension()) {
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
    const auto dimension = ct[0].dimension();
    const auto log_dimension = ct[0].log_dimension();
    const auto ct_mod_count = ct[0].component_count();
    const auto q_last = ct_moduli[ct_mod_count - 1]; // old last modulus
    const auto half_q_last = q_last / 2;
    const auto inv_t_mod_q_last = inverse_mod_prime(plain_modulus, q_last);
    std::vector<u64> q_last_mod_qi(ct_mod_count - 1);
    std::vector<u64> inv_q_last_mod_qi(ct_mod_count - 1);
    for (auto [q_i, q_last_reduced, inv_q_last] :
         zip(ct_moduli, q_last_mod_qi, inv_q_last_mod_qi)) {
        q_last_reduced = q_last % q_i; // need opt?
        inv_q_last = inverse_mod_prime(q_last, q_i);
    }

    for (auto &rns_poly : ct) {
        RnsPolynomial last_comp_copied{dimension, 1, std::vector{q_last}};
        last_comp_copied[0] = rns_poly[ct_mod_count - 1];
        intt_negacyclic_inplace_lazy(last_comp_copied);
        last_comp_copied *= inv_t_mod_q_last;
        batched_reduce_strict(q_last, dimension, last_comp_copied[0].data());
        // alias for clearness
        auto &last_comp_with_inv_t = last_comp_copied[0];

        RnsPolynomial subtract_part(dimension, ct_mod_count - 1, ct_moduli);
        for (auto [sub_part_comp, modulus, q_last_reduced] :
             zip(subtract_part, ct_moduli, q_last_mod_qi)) {
            // copy the last component and do reduction
            sub_part_comp = last_comp_with_inv_t;
            /* This reduction step needs further optimization. */
            batched_barrett(modulus, dimension, sub_part_comp.data());

            // The multiple of t needs to be of smallest possible abs value.
            for (auto [last_coeff_with_inv_t, sub_part_coeff] :
                 zip(last_comp_with_inv_t, sub_part_comp)) {
                if (last_coeff_with_inv_t >= half_q_last) {
                    sub_part_coeff += modulus - q_last_reduced;
                }
            }
        }
        subtract_part *= plain_modulus;
        ntt_negacyclic_inplace_lazy(subtract_part);

        rns_poly.remove_components();
        rns_poly -= subtract_part;
        rns_poly *= inv_q_last_mod_qi;
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
