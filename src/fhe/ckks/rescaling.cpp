#include "ckks.h"
#include "fhe/common/mod_arith.h"
#include "fhe/common/ntt.h"
#include "range/v3/view/zip.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>

using namespace ranges::views;

namespace hehub {
namespace ckks {

void rescale_by_one_prime_inplace(CkksCt &ct) {
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
    std::vector<u64> q_last_mod_qi(ct_mod_count - 1);
    std::vector<u64> inv_q_last_mod_qi(ct_mod_count - 1);
    for (auto [q_i, q_last_reduced, inv_q_last] :
         zip(ct_moduli, q_last_mod_qi, inv_q_last_mod_qi)) {
        q_last_reduced = q_last % q_i; // need opt?
        inv_q_last = inverse_mod_prime(q_last, q_i);
    }

    // this should be encapsulated as an RLWE utility in case useful in TFHE
    for (auto &rns_poly : ct) {
        RnsPolynomial last_comp{dimension, 1, std::vector{q_last}};
        last_comp[0] = rns_poly[ct_mod_count - 1];
        intt_negacyclic_inplace_lazy(last_comp);
        batched_reduce_strict(q_last, dimension, last_comp[0].data());

        auto &last_comp_coeffs = last_comp[0];
        RnsPolynomial remainder_q_last(dimension, ct_mod_count - 1, ct_moduli);
        for (auto [remainder_q_last_comp, modulus, q_last_reduced] :
             zip(remainder_q_last, ct_moduli, q_last_mod_qi)) {
            // copy the last component and do reduction
            remainder_q_last_comp = last_comp_coeffs;
            /* This reduction step needs further optimization. */
            batched_barrett(modulus, dimension, remainder_q_last_comp.data());

            // The remainder needs to be of smallest possible abs value, i.e. in
            // [-q_last/2, q_last/2).
            for (auto [last_comp_coeff, remainder_q_last_coeff] :
                 zip(last_comp_coeffs, remainder_q_last_comp)) {
                if (last_comp_coeff >= half_q_last) {
                    remainder_q_last_coeff += modulus - q_last_reduced;
                }
            }
        }
        ntt_negacyclic_inplace_lazy(remainder_q_last);

        rns_poly.remove_components();
        rns_poly -= remainder_q_last;
        rns_poly *= inv_q_last_mod_qi;
    }

    ct.scaling_factor /= q_last;
}

void rescale_inplace(CkksCt &ct, size_t dropping_primes) {
    if (dropping_primes == 1) {
        rescale_by_one_prime_inplace(ct);
    } else if (dropping_primes >= 2) {
        // TODO case
        throw "under development";
    } else {
        throw std::invalid_argument(
            "The number of primes to be dropped is not positive.");
    }
}

} // namespace ckks
} // namespace hehub
