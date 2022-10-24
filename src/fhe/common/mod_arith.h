/**
 * @file mod_arith.h
 * @brief Modular arithmetic.
 *
 */
#pragma once

#include "range/v3/view/zip.hpp"
#include "rns.h"
#include "type_defs.h"
#include <map>
#include <tuple>

namespace hehub {

void batched_barrett_lazy(const u64 modulus, const size_t vec_len, u64 vec[]);

inline void batched_barrett(const u64 modulus, const size_t vec_len,
                            u64 vec[]) {
    batched_barrett_lazy(modulus, vec_len, vec);

    for (size_t i = 0; i < vec_len; i++) {
        vec[i] -= (vec[i] >= modulus) ? modulus : 0;
    }
}

void batched_mul_mod_hybrid_lazy(const u64 modulus, const size_t vec_len,
                                 const u64 in_vec1[], const u64 in_vec2[],
                                 u64 out_vec[]);

inline void batched_mul_mod_hybrid(const u64 modulus, const size_t vec_len,
                                   const u64 in_vec1[], const u64 in_vec2[],
                                   u64 out_vec[]) {
    batched_mul_mod_hybrid_lazy(modulus, vec_len, in_vec1, in_vec2, out_vec);

    for (size_t i = 0; i < vec_len; i++) {
        out_vec[i] -= (out_vec[i] >= modulus) ? modulus : 0;
    }
}

void batched_mul_mod_barrett_lazy(const u64 modulus, const size_t vec_len,
                                  const u64 in_vec1[], const u64 in_vec2[],
                                  u64 out_vec[]);

inline void batched_mul_mod_barrett(const u64 modulus, const size_t vec_len,
                                    const u64 in_vec1[], const u64 in_vec2[],
                                    u64 out_vec[]) {
    batched_mul_mod_barrett_lazy(modulus, vec_len, in_vec1, in_vec2, out_vec);

    for (size_t i = 0; i < vec_len; i++) {
        out_vec[i] -= (out_vec[i] >= modulus) ? modulus : 0;
    }
}

void batched_montgomery_128_lazy(const u64 modulus, const size_t len,
                                 const u128 in[], u64 out[]);

inline void batched_reduce_strict(const u64 modulus, const size_t vec_len,
                                  u64 vec[]) {
    for (size_t i = 0; i < vec_len; i++) {
        vec[i] -= (vec[i] >= modulus) ? modulus : 0;
    }
}

inline void reduce_strict(RnsPolynomial &rns_poly) {
    const auto &moduli = rns_poly.modulus_vec();
    const auto dimension = rns_poly.dimension();

    for (auto [component, modulus] : ranges::views::zip(rns_poly, moduli)) {
        batched_reduce_strict(modulus, dimension, component.data());
    }
}

inline u64 mul_mod_harvey_lazy(const u64 modulus, const u64 in1, const u64 in2,
                               const u64 in2_harvey) {
    u64 approx_quotient = (u128)in1 * in2_harvey >> 64;
    return (u128)in1 * in2 - (u128)approx_quotient * modulus;
}

extern std::map<std::pair<u64, u64>, u64> modular_inverse_table;

u64 inverse_mod_prime(const u64 elem, const u64 prime);

} // namespace hehub
