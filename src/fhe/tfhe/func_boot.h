/**
 * @file func_boot.h
 * @brief TODO
 *
 */

#pragma once

#include "fhe/common/rns.h"
#include "fhe/primitives/lwe.h"
#include "fhe/primitives/rgsw.h"

namespace hehub {
namespace tfhe {

/**
 * @brief Perform a functional bootstrapping on an input LWE ciphertext with a
 * lookup table (LUT) in form of a polynomial, and the bootstrapping keys.
 * @note The mechanism of functional bootstrapping (FBS) is to iteratively
 * compute the polynomial lut_poly * X^(b + a_0*s_0 + a_1*s_1 + ...) under
 * encryption, and then extract an LWE ciphertext carrying the required message
 * from an RLWE ciphertext. The iteration is also refered to as blind rotation,
 * which makes blind selections each round based on the value of s_i. After the
 * blind rotation, a specific coefficient of lut_poly is move to the constant
 * term of the resulting polynomial. This resulting polynomial is under
 * encryption of RLWE, from which an LWE ciphertext that carries the constant
 * term can be extracted.
 * @note The functional bootstrapping procedure is negacyclic, which means the
 * resulting message can be opposite to what is encoded in the coefficients of
 * lut_poly. The procedure of fully functional bootstrapping below avoids this.
 * @note For more details of functional bootstrapping, see Chillotti et al.
 * (https://doi.org/10.1007/s00145-019-09319-x).
 * @param ct The input LWE ciphertext.
 * @param lut_poly A lookup table in form of a polynomial, which encodes the
 * function to evaluate.
 * @param bootstrap_keys The bootstrapping keys, i.e. the RGSW encryption on
 * entries of the secret key.
 * @return LweCt
 */
LweCt functional_bootstrap(const LweCt &ct, const RnsPolynomial &lut_poly,
                           const std::vector<RgswCt> &bootstrap_keys);

/**
 * @brief Get the redundant most-significant-bit of the message of ct which is
 * supposed to go through a functional bootstrapping next.
 * @note Let n be the dimension of ct. When doing a blind rotation, it is
 * required ct's modulus is n as well. Now notice the blind rotation procedure
 * does almost an LWE decryption on the exponent of X, just except that the
 * numbers on the exponent are modulo 2n. Let ct = (a,b) and the secret key = s,
 * and let <a,s>+b = m+e (mod n), then <a,s>+b = m+e+kn (mod 2n) where k can be
 * 0 or 1. The k is refered to as the redundant most-significant-bit, and due to
 * the negacyclicity property of polynomials it will mess the result if it is
 * not 0.
 * @note This function runs a FBS on ct, but with a specifically designed LUT
 * polynomial, so that it can utilize the negacyclicity property for obtaining
 * k, cf. Yang et al. (https://eprint.iacr.org/2021/1347.pdf) or Liu et al.
 * (https://eprint.iacr.org/2021/1337.pdf).
 * @param ct The input LWE ciphertext.
 * @param bootstrap_keys The bootstrapping keys, i.e. the RGSW encryption on
 * entries of the secret key.
 * @return LweCt
 */
LweCt get_redundant_msb(const LweCt &ct,
                        const std::vector<RgswCt> &bootstrap_keys);

/**
 * @brief Perform a fully functional bootstrapping on an input LWE ciphertext
 * with a lookup table (LUT) in form of a polynomial, and the bootstrapping
 * keys.
 * @note The procedure of fully functional bootstrapping (FFBS) avoids the
 * negacyclicity problem of FBS. It first invokes get_redundant_msb on ct,
 * obtaining an LWE ciphertext of kn (when modulus is 2n), where k is the
 * redundant MSB. Then subtracting this ciphertext from ct ensures that we
 * obtain a ciphertext without the redundant MSB when performing blind rotation,
 * on which we can finally run FBS with the LUT polynomial.
 * @note The cost of efficiency is about 2x compared to a FBS, since FFBS
 * basically makes two invocation of FBS, one for getting redundant MSB and the
 * other for looking up required message from lut_poly. cf. Yang et al.
 * (https://eprint.iacr.org/2021/1347.pdf) or Liu et al.
 * (https://eprint.iacr.org/2021/1337.pdf).
 * @param ct The input LWE ciphertext.
 * @param lut_poly A lookup table in form of a polynomial, which encodes the
 * function to evaluate.
 * @param bootstrap_keys The bootstrapping keys, i.e. the RGSW encryption on
 * entries of the secret key.
 * @return LweCt
 */
LweCt fully_functional_bootstrap(const LweCt &ct, const RnsPolynomial &lut_poly,
                                 const std::vector<RgswCt> &bootstrap_keys);

} // namespace tfhe
} // namespace hehub
