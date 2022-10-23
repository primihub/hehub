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

/**
 * @brief TODO
 *
 * @param modulus
 * @param vec_len
 * @param vec
 */
void batched_barrett_lazy(const u64 modulus, const size_t vec_len, u64 vec[]);

/**
 * @brief TODO
 *
 * @param modulus
 * @param vec_len
 * @param vec
 */
inline void batched_barrett(const u64 modulus, const size_t vec_len,
                            u64 vec[]) {
    batched_barrett_lazy(modulus, vec_len, vec);

    for (size_t i = 0; i < vec_len; i++) {
        vec[i] -= (vec[i] >= modulus) ? modulus : 0;
    }
}

/**
 * @brief An optimized method for calculating modular multiplication (between
 * vectors for efficiency) for a modulus of less than 64-bit, and input vectors
 * arbituary (non-fixed). Here we apply a hybrid method, which does a montgomery
 * reduction on in_vec1[i] * in_vec2[i] producing a result equal to (in_vec1[i]
 * * in_vec2[i] * C) % modulus where C is a montgomery constant related to the
 * modulus, and then does a harvey reduction to remove C. (cf.
 * https://en.wikipedia.org/wiki/Montgomery_modular_multiplication and
 * https://doi.org/10.1016/j.jsc.2013.09.002) The resulting out_vec[i]'s are
 * left in [0, 2*modulus).
 *
 * @param[in] modulus The modulus.
 * @param[in] vec_len The length of input vectors.
 * @param[in] in_vec1 Input vector of length vec_len
 * @param[in] in_vec2 Input vector of length vec_len
 * @param[out] out_vec Output vector s.t. out_vec[i] == in_vec1[i] * in_vec2[i]
 * % modulus.
 */
void batched_mul_mod_hybrid_lazy(const u64 modulus, const size_t vec_len,
                                 const u64 in_vec1[], const u64 in_vec2[],
                                 u64 out_vec[]);

/**
 * @brief An optimized method for calculating modular multiplication (between
 * vectors for efficiency) for a modulus of less than 64-bit, and input vectors
 * arbituary (non-fixed). Here we apply a hybrid method, which does a montgomery
 * reduction on in_vec1[i] * in_vec2[i] producing a result equal to (in_vec1[i]
 * * in_vec2[i] * C) % modulus where C is a montgomery constant related to the
 * modulus, and then does a harvey reduction to remove C. (cf.
 * https://en.wikipedia.org/wiki/Montgomery_modular_multiplication and
 * https://doi.org/10.1016/j.jsc.2013.09.002)
 *
 * @param[in] modulus The modulus.
 * @param[in] vec_len The length of input vectors.
 * @param[in] in_vec1 Input vector of length vec_len
 * @param[in] in_vec2 Input vector of length vec_len
 * @param[out] out_vec Output vector s.t. out_vec[i] == in_vec1[i] * in_vec2[i]
 * % modulus.
 */
inline void batched_mul_mod_hybrid(const u64 modulus, const size_t vec_len,
                                   const u64 in_vec1[], const u64 in_vec2[],
                                   u64 out_vec[]) {
    batched_mul_mod_hybrid_lazy(modulus, vec_len, in_vec1, in_vec2, out_vec);

    for (size_t i = 0; i < vec_len; i++) {
        out_vec[i] -= (out_vec[i] >= modulus) ? modulus : 0;
    }
}

/**
 * @brief Use barrett reduction for calculating modular multiplication (between
 * vectors for efficiency) for a modulus of less than 64-bit, and input vectors
 * arbituary (non-fixed). After the product of in_vec1[i] and in_vec2[i] is
 * produced, barrett reduction works as follows: A barrett constant c is
 * pre-computed as floor(2**128 / modulus). Let the product of in_vec1[i] and
 * in_vec2[i] is a, then floor(a / modulus) ≈ floor(a * c / 2**128), from which
 * we can obtain a % modulus with possibly one additional multiple of modulus.
 * The resulting out_vec[i]'s are left in [0, 2*modulus).
 * @see https://en.wikipedia.org/wiki/Barrett_reduction
 *
 * @param[in] modulus The modulus.
 * @param[in] vec_len The length of input vectors.
 * @param[in] in_vec1 Input vector of length vec_len
 * @param[in] in_vec2 Input vector of length vec_len
 * @param[out] out_vec Output vector s.t. out_vec[i] == in_vec1[i] * in_vec2[i]
 * % modulus.
 */
void batched_mul_mod_barrett_lazy(const u64 modulus, const size_t vec_len,
                                  const u64 in_vec1[], const u64 in_vec2[],
                                  u64 out_vec[]);

/**
 * @brief Use barrett reduction for calculating modular multiplication (between
 * vectors for efficiency) for a modulus of less than 64-bit, and input vectors
 * arbituary (non-fixed). After the product of in_vec1[i] and in_vec2[i] is
 * produced, barrett reduction works as follows: A barrett constant c is
 * pre-computed as floor(2**128 / modulus). Let the product of in_vec1[i] and
 * in_vec2[i] is a, then floor(a / modulus) ≈ floor(a * c / 2**128), from which
 * we can obtain a % modulus with possibly one additional multiple of modulus,
 * cf. https://en.wikipedia.org/wiki/Barrett_reduction.
 *
 * @param[in] modulus The modulus.
 * @param[in] vec_len The length of input vectors.
 * @param[in] in_vec1 Input vector of length vec_len
 * @param[in] in_vec2 Input vector of length vec_len
 * @param[out] out_vec Output vector s.t. out_vec[i] == in_vec1[i] * in_vec2[i]
 * % modulus.
 */
inline void batched_mul_mod_barrett(const u64 modulus, const size_t vec_len,
                                    const u64 in_vec1[], const u64 in_vec2[],
                                    u64 out_vec[]) {
    batched_mul_mod_barrett_lazy(modulus, vec_len, in_vec1, in_vec2, out_vec);

    for (size_t i = 0; i < vec_len; i++) {
        out_vec[i] -= (out_vec[i] >= modulus) ? modulus : 0;
    }
}

/**
 * @brief TODO
 *
 * @param modulus
 * @param len
 * @param in
 * @param out
 */
void batched_montgomery_128_lazy(const u64 modulus, const size_t len,
                                 const u128 in[], u64 out[]);

/**
 * @brief TODO
 *
 * @param modulus
 * @param vec_len
 * @param vec
 */
inline void batched_reduce_strict(const u64 modulus, const size_t vec_len,
                                  u64 vec[]) {
    for (size_t i = 0; i < vec_len; i++) {
        vec[i] -= (vec[i] >= modulus) ? modulus : 0;
    }
}

/**
 * @brief TODO
 *
 * @param rns_poly
 */
inline void reduce_strict(RnsPolynomial &rns_poly) {
    const auto &moduli = rns_poly.modulus_vec();
    const auto dimension = rns_poly.dimension();

    for (auto [component, modulus] : ranges::views::zip(rns_poly, moduli)) {
        batched_reduce_strict(modulus, dimension, component.data());
    }
}

/**
 * @brief TODO
 *
 * @param modulus
 * @param in1
 * @param in2
 * @param in2_harvey
 * @return u64
 */
inline u64 mul_mod_harvey_lazy(const u64 modulus, const u64 in1, const u64 in2,
                               const u64 in2_harvey) {
    u64 approx_quotient = (u128)in1 * in2_harvey >> 64;
    return (u128)in1 * in2 - (u128)approx_quotient * modulus;
}

/**
 * @brief TODO
 *
 */
extern std::map<std::pair<u64, u64>, u64> modular_inverse_table;

/**
 * @brief
 *
 * @param elem
 * @param prime
 * @return u64
 */
u64 inverse_mod_prime(const u64 elem, const u64 prime);

} // namespace hehub
