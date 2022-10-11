/**
 * @file permutation.h
 * @brief TODO
 *
 */

#pragma once
#include "rnspolynomial.h"
#include "type_defs.h"

namespace hehub {

/**
 * @brief TODO
 *
 * @param x
 * @param bit_len
 * @return u64
 */
inline u64 __bit_rev_naive(u64 x, int bit_len) {
#ifdef HEHUB_DEBUG
    if (bit_len < 0 || bit_len > 64) {
        throw std::invalid_argument("bit_len");
    }
    if (bit_len != 64 && x >= (1ULL << bit_len)) {
        throw std::invalid_argument("x");
    }
#endif
    if (bit_len <= 1) {
        return x;
    }
    u64 mh = (1 << (bit_len - 1));
    u64 ml = 1;
    for (; mh >= ml; mh >>= 1, ml <<= 1) {
        if ((mh & x) && (ml & x)) {
        } else if (!(mh & x) && !(ml & x)) {
        } else if (!(mh & x) && (ml & x)) {
            x |= mh;
            x ^= ml;
        } else {
            x |= ml;
            x ^= mh;
        }
    }
    return x;
}

/**
 * @brief TODO
 *
 * @param x
 * @param bit_len
 * @return u64
 */
inline u64 __bit_rev_naive_16(u64 x, int bit_len) {
#ifdef HEHUB_DEBUG
    if (bit_len < 0 || bit_len > 16) {
        throw std::invalid_argument("bit_len");
    }
    if (x >= (1ULL << bit_len)) {
        throw std::invalid_argument("x");
    }
#endif
    x = ((x & 0xFF00FF00) >> 8) | ((x & 0x00FF00FF) << 8);
    x = ((x & 0xF0F0F0F0) >> 4) | ((x & 0x0F0F0F0F) << 4);
    x = ((x & 0xCCCCCCCC) >> 2) | ((x & 0x33333333) << 2);
    x = ((x & 0xAAAAAAAA) >> 1) | ((x & 0x55555555) << 1);
    return x >> (16 - bit_len);
}

/// @brief The exponents of the Galois transformations when interpreted as
/// exponentiation operations.
/// @note The Galois group of Q(ξ)/Q with ξ being a 2-power m-th primitive root
/// is of structure Z_{m/4}⊕Z_2 when m >= 8, where the elements (the Galois
/// transformations) are exponentiations with exponents 1, 3, 5, ... and m-1
/// respectively, some of which are of maximum possible order m/4. In the
/// implementation we choose such a generator 3 since it is small and universal,
/// i.e. fits any m >= 8 (which can be demonstrated by induction).
/// @note Since the negecyclic Fourier transformations produce values of a
/// polynomial on unity roots ξ, ξ^3, ..., ξ^(m-1) with ξ being a primitive
/// root, an exponentiation modifies the index of the corresponding roots resp.
/// to ξ by multiplying them with a factor. Hence here we refer to the exponents
/// as root index factors.
/// @note This table can also be viewed as the data encoding order when using
/// SIMD encoding (after minor processing such as bit reversal), since a next
/// datum should be place in the rotated position of the previous one.
std::vector<u32> &root_index_factors();

/// @brief TODO
/// @param loglen
/// @return
const std::vector<size_t> &dlog_mod_2power_table(size_t loglen);

/**
 * @brief Perform a cycle on the values of the polynomial, which is induced from
 * the Galois transformation from a cyclic subgroup of order poly_len / 2.
 * @note The Galois group of Q(ξ)/Q with ξ being a 2-power m-th primitive unity
 * root is of structure Z_{m/4}⊕Z_2, hence a Galois transformations of Q(ξ)/Q is
 * either a cycle from a subgroup of order m/4, or an involution. On the other
 * hand, when q is an unramified prime in Q(ξ) (or product of such primes), the
 * ring R_q = Z_q[X]/(X^n + 1) is isomorphic to an n-dimensional q-linear space,
 * whose automorphism group is the group of all n-dimensional q-linear
 * transformations. It is homomorphic encryption that requires the automorphism
 * to be bounded (resp. to inf. norm). Hence we cannot utilize arbituary
 * q-linear transformations but only those induced by the Galois ones.
 * @param poly_ntt An RnsPolynomial in NTT value form.
 * @param step The cycle step resp. to the (poly_len / 2)-ordered subgroup.
 * @return RnsPolynomial
 */
RnsPolynomial cycle(const RnsPolynomial &poly_ntt, const size_t step);

/**
 * @brief Perform an involution on the values of the polynomial, which is
 * induced from the Galois transformation of conjugation.
 * @note The Galois group of Q(ξ)/Q with ξ being a 2-power m-th primitive unity
 * root is of structure Z_{m/4}⊕Z_2, hence a Galois transformations of Q(ξ)/Q is
 * either a cycle from a subgroup of order m/4, or an involution. On the other
 * hand, when q is an unramified prime in Q(ξ) (or product of such primes), the
 * ring R_q = Z_q[X]/(X^n + 1) is isomorphic to an n-dimensional q-linear space,
 * whose automorphism group is the group of all n-dimensional q-linear
 * transformations. It is homomorphic encryption that requires the automorphism
 * to be bounded (resp. to inf. norm). Hence we cannot utilize arbituary
 * q-linear transformations but only those induced by the Galois ones.
 * @note There are three involutions in the Galois transformation set. This
 * function performs the one that corresponds to complex conjugation.
 * @param poly_ntt An RnsPolynomial in NTT value form.
 * @return RnsPolynomial
 */
RnsPolynomial involute(const RnsPolynomial &poly_ntt);

} // namespace hehub
