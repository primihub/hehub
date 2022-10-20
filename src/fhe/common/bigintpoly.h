#pragma once

#include "bigint.h"
#include "rns.h"

namespace hehub {

/**
 * @class UBigIntPoly
 * @brief Polynomial with coefficients implemented as big integers.
 */
class UBigIntPoly {
public:
    /**
     * @brief Construct a new UBigIntPoly object by converting from an RNS
     * polynomial.
     * @param rns_poly The source RNS polynomial.
     */
    UBigIntPoly(const RnsPolynomial &rns_poly);

    /**
     * @brief Get the length of this polynomial, which is its degree + 1.
     * @return const size_t 
     */
    inline const size_t dimension() const { return coeffs_.size(); }

    /**
     * @brief Get the i-th coefficient of this polynomial.
     * @param i The specified coefficient index.
     * @return UBInt & 
     */
    inline UBInt &operator[](int i) { return coeffs_[i]; }

    /**
     * @brief Get the i-th coefficient of this polynomial as read-only.
     * @param i The specified coefficient index.
     * @return const UBInt & 
     */
    inline const UBInt &operator[](int i) const { return coeffs_[i]; }

    /**
     * @brief Print the unsigned big-int polynomial.
     * @return std::ostream & 
     */
    friend std::ostream &operator<<(std::ostream &, const UBigIntPoly &);

private:
    /// A vector storing all the coefficients.
    std::vector<UBInt> coeffs_;
};

} // namespace hehub
