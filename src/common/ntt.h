/**
 * @file ntt.h
 * @brief Number Theoretic Transforms that transform a polynomial modulo both of
 * a prime q and (X^n + 1) for n being a power of 2 (a.k.a being an element of
 * the so-called nega-cyclic ring Z_q[X]/(X^n + 1)) between coefficient and
 * value form of presentation, cf.
 * https://en.wikipedia.org/wiki/Discrete_Fourier_transform_over_a_ring.
 */

#pragma once

#include "type_defs.h"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace hehub {

/**
 * @brief The function carries out forward NTT operation inplace, the input and
 * output of which is an element of the negacylic ring Z_q[X]/(X^n + 1) where q
 * and log2(n) are provided from parameters. This element is viewed as a
 * polynomial in coefficient form. The output values are left in [0, 2*modulus).
 * @param[in] modulus The modulus q.
 * @param[in] log_poly_len The log value of the length of the polynomial n, i.e.
 * the number of coefficients.
 * @param[inout] coeffs The polynomial in coefficient form.
 */
void ntt_negacyclic_inplace_lazy(const u64 modulus, const size_t log_poly_len,
                                 u64 coeffs[]);

/**
 * @brief The function carries out inverse NTT operation inplace, the input and
 * output of which is an element of the negacylic ring Z_q[X]/(X^n + 1) where q
 * and log2(n) are provided from parameters. This element is viewed as the
 * values evaluated from itself viewed as an integral polynomial modulo q.  The
 * output values are left in [0, 2*modulus).
 * @param[in] modulus The modulus q.
 * @param[in] log_poly_len The log value of the length of the polynomial n,
 * which is also the number of values.
 * @param[inout] values The values evaluated from the polynomial.
 */
void intt_negacyclic_inplace_lazy(const u64 modulus, const size_t log_poly_len,
                                  u64 values[]);

/**
 * @brief Create the NTT roots for all the input moduli immediately, which will
 * be used in the process of NTT and INTT. N.B., this functionality is
 * _optional_ since we use a lazy strategy in creating these NTT roots, i.e.
 * create and store them when being used the first time.
 */
void cache_ntt_factors_strict(const std::vector<u64> &moduli,
                              const u64 log_poly_len);

} // namespace hehub
