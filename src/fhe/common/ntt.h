/**
 * @file ntt.h
 * @brief Number Theoretic Transforms that transform a polynomial modulo both of
 * a prime q and (X^n + 1) for n being a power of 2 (a.k.a being an element of
 * the negacyclic ring Z_q[X]/(X^n + 1)) between coefficient and value form of
 * presentation.
 * @note cf.
 * https://en.wikipedia.org/wiki/Discrete_Fourier_transform_over_a_ring.
 */

#pragma once

#include "mod_arith.h"
#include "range/v3/view/zip.hpp"
#include "rns.h"
#include "type_defs.h"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace hehub {

/**
 * @brief The function carries out forward NTT operation inplace, the input and
 * output of which is an element of the negacyclic ring Z_q[X]/(X^n + 1) where q
 * and log2(n) are provided from parameters. This element is viewed as a
 * polynomial in coefficient form. The output values are left in [0, 2*modulus).
 * @param[in] log_dimension The log value of the length of the polynomial, i.e.
 * the number of coefficients.
 * @param[in] modulus The modulus q.
 * @param[inout] coeffs The polynomial in coefficient form.
 */
void ntt_negacyclic_inplace_lazy(const size_t log_dimension, const u64 modulus,
                                 u64 coeffs[]);

/**
 * @brief TODO
 *
 * @param[inout] rns_poly
 */
inline void ntt_negacyclic_inplace_lazy(RnsPolynomial &rns_poly) {
    const auto component_count = rns_poly.component_count();
    const auto log_dimension = rns_poly.log_dimension();
    const auto &moduli = rns_poly.modulus_vec();

    for (auto [component, modulus] : ranges::views::zip(rns_poly, moduli)) {
        ntt_negacyclic_inplace_lazy(log_dimension, modulus, component.data());
    }

    rns_poly.rep_form = PolyRepForm::value;
}

/**
 * @brief The function carries out inverse NTT operation inplace, the input and
 * output of which is an element of the negacyclic ring Z_q[X]/(X^n + 1) where q
 * and log2(n) are provided from parameters. This element is viewed as the
 * values evaluated from itself viewed as an integral polynomial modulo q.  The
 * output values are left in [0, 2*modulus).
 * @param[in] log_dimension The log value of the length of the polynomial,
 * which is the number of input values.
 * @param[in] modulus The modulus q.
 * @param[inout] values The values evaluated from the polynomial.
 */
void intt_negacyclic_inplace_lazy(const size_t log_dimension, const u64 modulus,
                                  u64 values[]);

/**
 * @brief TODO
 *
 * @param[inout] rns_poly
 */
inline void intt_negacyclic_inplace_lazy(RnsPolynomial &rns_poly) {
    const auto component_count = rns_poly.component_count();
    const auto log_dimension = rns_poly.log_dimension();
    const auto &moduli = rns_poly.modulus_vec();

    for (auto [component, modulus] : ranges::views::zip(rns_poly, moduli)) {
        intt_negacyclic_inplace_lazy(log_dimension, modulus, component.data());
    }

    rns_poly.rep_form = PolyRepForm::coeff;
}

/**
 * @brief TODO
 *
 * @param[inout] rns_poly
 */
inline void intt_negacyclic_inplace(RnsPolynomial &rns_poly) {
    intt_negacyclic_inplace_lazy(rns_poly);
    reduce_strict(rns_poly);
}

/**
 * @brief Create the NTT factors for all the input moduli immediately, which
 * will be used in the process of NTT and INTT. N.B., this functionality is
 * _optional_ in use since we use a lazy strategy in creating these NTT factors,
 * i.e. create and store them when needed the first time.
 * @param[in] log_dimension The log value of the length of the polynomial.
 * @param[in] moduli The moduli modulo which the NTT factors are produced.
 */
void cache_ntt_factors_strict(const u64 log_dimension,
                              const std::vector<u64> &moduli);

} // namespace hehub
