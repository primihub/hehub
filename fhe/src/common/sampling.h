/**
 * @file sampling.h
 * @brief Generate random polynomials.
 */

#include "rnspolynomial.h"

namespace hehub {

/**
 * @brief Get a random ternary RnsPolynomial object with specified dimensions,
 * which represents a polynomial with coefficients uniform from {-1, 0, 1}. The
 * result will be in NTT form.
 * @param poly_dim Polynomial dimensions for initializing the RnsPolynomial.
 * @return RnsPolynomial
 */
RnsPolynomial get_rand_ternary_poly(const PolyDimensions &poly_dim);

/**
 * @brief Get a random uniform RnsPolynomial object with specified dimensions,
 * where the k-th component polynomial's coefficients will be uniform from {0,
 * ..., moduli[k] - 1}.
 * @param poly_dim Polynomial dimensions for initializing the RnsPolynomial.
 * @param form Required representation form of the resulting polynomial.
 * @return RnsPolynomial
 */
RnsPolynomial get_rand_uniform_poly(const PolyDimensions &poly_dim,
                                    PolyRepForm form = PolyRepForm::coeff);

/**
 * @brief Get a random gaussian RnsPolynomial object with specified dimensions,
 * which represents a polynomial with coefficients sampled (and rounded) from
 * a gaussian distribution. The result will be in NTT form.
 * @param poly_dim Polynomial dimensions for initializing the RnsPolynomial.
 * @param std_dev Standard deviation of the gaussian distribution.
 * @return RnsPolynomial
 */
RnsPolynomial get_rand_gaussian_poly(const PolyDimensions &poly_dim,
                                     double std_dev = 3.2);

/**
 * @brief TODO
 *
 * @param poly_dim
 * @param form
 * @return RnsPolynomial
 */
RnsPolynomial get_zero_poly(const PolyDimensions &poly_dim,
                            PolyRepForm form = PolyRepForm::value);

} // namespace hehub
