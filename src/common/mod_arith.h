/**
 * @file mod_arith.h
 * @brief Modular arithmetic.
 *
 */
#pragma once

#include "rnspolynomial.h"
#include "type_defs.h"
#include <tuple>

namespace hehub {

using MulModLUT = std::tuple<u64, u64, u64>;

/**
 * @brief An optimized method for calculating modular multiplication (between
 * vectors for efficiency) for a modulus of less than 64-bit, and input vectors
 * arbituary (non-fixed). Here we apply a hybrid method, which does a montgomery
 * reduction on f[i] * g[i] producing a result equal to (f[i] * g[i] * C) %
 * modulus where C is a montgomery constant related to the modulus, and then
 * does a harvey reduction to remove C. (cf.
 * https://en.wikipedia.org/wiki/Montgomery_modular_multiplication and
 * https://doi.org/10.1016/j.jsc.2013.09.002)
 *
 * @param[in] modulus The modulus.
 * @param[in] vec_len The length of input vectors.
 * @param[in] f Input vector of length vec_len
 * @param[in] g Input vector of length vec_len
 * @param[out] h Output vector s.t. h[i] == f[i] * g[i] % modulus.
 */
void vector_mul_mod_hybrid(const u64 modulus, const size_t vec_len, const u64 *f,
                    const u64 *g, u64 *h);

/**
 * @brief An optimized method for calculating modular multiplication (between
 * vectors for efficiency) for a modulus of less than 64-bit, and input vectors
 * arbituary (non-fixed). Here we apply a hybrid method, which does a montgomery
 * reduction on f[i] * g[i] producing a result equal to (f[i] * g[i] * C) %
 * modulus where C is a montgomery constant related to the modulus, and then
 * does a harvey reduction to remove C. (cf.
 * https://en.wikipedia.org/wiki/Montgomery_modular_multiplication and
 * https://doi.org/10.1016/j.jsc.2013.09.002)
 *
 * @param[in] modulus The modulus.
 * @param[in] vec_len The length of input vectors.
 * @param[in] f Input vector of length vec_len
 * @param[in] g Input vector of length vec_len
 * @param[out] h Output vector s.t. h[i] == f[i] * g[i] % modulus.
 */
inline void vector_mul_mod_hybrid(const u64 modulus, const size_t vec_len,
                           const RnsPolynomial::ComponentData &f,
                           const RnsPolynomial::ComponentData &g,
                           RnsPolynomial::ComponentData &h) {
    vector_mul_mod_hybrid(modulus, vec_len, f.data(), g.data(), h.data());
}

/**
 * @brief Use barrett reduction for calculating modular multiplication (between
 * vectors for efficiency) for a modulus of less than 64-bit, and input vectors
 * arbituary (non-fixed). After the product of f[i] and g[i] is produced,
 * barrett reduction works as follows: A barrett constant c is pre-computed as
 * floor(2**128 / modulus). Let the product of f[i] and g[i] is a, then floor(a
 * / modulus) ≈ floor(a * c / 2**128), from which we can obtain a % modulus with
 * possibly one additional multiple of modulus, cf.
 * https://en.wikipedia.org/wiki/Barrett_reduction.
 *
 * @param[in] modulus The modulus.
 * @param[in] vec_len The length of input vectors.
 * @param[in] f Input vector of length vec_len
 * @param[in] g Input vector of length vec_len
 * @param[out] h Output vector s.t. h[i] == f[i] * g[i] % modulus.
 */
void vector_mul_mod_barrett(const u64 modulus, const size_t vec_len,
                            const u64 *f, const u64 *g, u64 *h);
                            
/**
 * @brief Use barrett reduction for calculating modular multiplication (between
 * vectors for efficiency) for a modulus of less than 64-bit, and input vectors
 * arbituary (non-fixed). After the product of f[i] and g[i] is produced,
 * barrett reduction works as follows: A barrett constant c is pre-computed as
 * floor(2**128 / modulus). Let the product of f[i] and g[i] is a, then floor(a
 * / modulus) ≈ floor(a * c / 2**128), from which we can obtain a % modulus with
 * possibly one additional multiple of modulus, cf.
 * https://en.wikipedia.org/wiki/Barrett_reduction.
 *
 * @param[in] modulus The modulus.
 * @param[in] vec_len The length of input vectors.
 * @param[in] f Input vector of length vec_len
 * @param[in] g Input vector of length vec_len
 * @param[out] h Output vector s.t. h[i] == f[i] * g[i] % modulus.
 */
inline void vector_mul_mod_barrett(const u64 modulus, const size_t vec_len,
                           const RnsPolynomial::ComponentData &f,
                           const RnsPolynomial::ComponentData &g,
                           RnsPolynomial::ComponentData &h) {
    vector_mul_mod_barrett(modulus, vec_len, f.data(), g.data(), h.data());
}

} // namespace hehub
