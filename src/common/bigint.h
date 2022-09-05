/**
 * @file ubigint.h
 * @brief Class of big integers and utility of CRT composing.
 *
 */

#pragma once

#include "type_defs.h"
#include <cstring>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace hehub {

/**
 * @class UBInt
 * @brief Class of unsigned big integers which supports basic arithmetic,
 * comparison and printing.
 *
 */
class UBInt {
public:
    /**
     * @brief Construct a new unsigned big integer object from a unsigned
     * integer.
     * @param nr The source integer.
     */
    UBInt(u64 nr = 0);

    /**
     * @brief Construct a new unsigned big integer object from a string
     * specifying the digits.
     * @param str The source string.
     */
    UBInt(const std::string &str);

    /**
     * @brief Construct a new unsigned big integer object from a C-string
     * specifying the digits.
     * @param str The source C-string.
     */
    UBInt(const char *str);

    /**
     * @brief Construct a new unsigned big integer object by copying another
     * one.
     * @param other The source UBInt object.
     */
    UBInt(const UBInt &other) = default;

    /**
     * @brief Construct a new unsigned big integer object by moving in another
     * one.
     * @param other The source UBInt object.
     */
    UBInt(UBInt &&other) = default;

    /**
     * @brief Direct assignment by copying another one.
     * @param other
     * @return UBInt &
     */
    UBInt &operator=(const UBInt &other) = default;

    /**
     * @brief Direct assignment by moving in another one.
     * @param other
     * @return UBInt &
     */
    UBInt &operator=(UBInt &&other) = default;

    /**
     * @brief Pre-incrementation.
     */
    UBInt &operator++();

    /**
     * @brief Post-incrementation.
     */
    UBInt operator++(int _dummy_);

    /**
     * @brief Pre-decrementation.
     */
    UBInt &operator--();

    /**
     * @brief Post-decrementation.
     */
    UBInt operator--(int _dummy_);

    /**
     * @brief Add b to a, and return the result.
     * @param a
     * @param b
     * @return UBInt &
     */
    friend UBInt &operator+=(UBInt &a, const UBInt &b);

    /**
     * @brief Return the sum of a and b.
     * @param a
     * @param b
     * @return UBInt
     */
    friend UBInt operator+(const UBInt &a, const UBInt &b);

    /**
     * @brief Subtract a by b, and return the result.
     * @param a
     * @param b
     * @return UBInt &
     */
    friend UBInt &operator-=(UBInt &a, const UBInt &b);

    /**
     * @brief Return the difference between a and b.
     * @param a
     * @param b
     * @return UBInt
     */
    friend UBInt operator-(const UBInt &a, const UBInt &b);

    /**
     * @brief Multiply a by b, and return the result.
     * @param a
     * @param b
     * @return UBInt&
     */
    friend UBInt &operator*=(UBInt &a, const UBInt &b);

    /**
     * @brief Return the product of a and b.
     * @param a
     * @param b
     * @return UBInt
     */
    friend UBInt operator*(const UBInt &a, const UBInt &b);

    /**
     * @brief Divide a by b, and return the result.
     * @param a
     * @param b
     * @return UBInt &
     */
    friend UBInt &operator/=(UBInt &a, const UBInt &b);

    /**
     * @brief Return the result of a divided by b.
     * @param a
     * @param b
     * @return UBInt
     */
    friend UBInt operator/(const UBInt &a, const UBInt &b);

    /**
     * @brief Reduce a modulo b, and return the result.
     * @param a
     * @param b
     * @return UBInt &
     */
    friend UBInt &operator%=(UBInt &a, const UBInt &b);

    /**
     * @brief Return the result of a reduced modulo b.
     * @param a
     * @param b
     * @return UBInt
     */
    friend UBInt operator%(const UBInt &a, const UBInt &b);

    /**
     * @brief Square root function.
     * @return UBInt
     */
    friend UBInt sqrt(const UBInt &);

    /**
     * @brief Tell if a is equivalent to b.
     * @param a
     * @param b
     * @return bool
     */
    friend bool operator==(const UBInt &a, const UBInt &b);

    /**
     * @brief Tell if a is non-equivalent to b.
     * @param a
     * @param b
     * @return bool
     */
    friend bool operator!=(const UBInt &a, const UBInt &b);

    /**
     * @brief Tell if a is greater than b.
     * @param a
     * @param b
     * @return bool
     */
    friend bool operator>(const UBInt &a, const UBInt &b);

    /**
     * @brief Tell if a is greater than or equal to b.
     * @param a
     * @param b
     * @return bool
     */
    friend bool operator>=(const UBInt &a, const UBInt &b);

    /**
     * @brief Tell if a is less than b.
     * @param a
     * @param b
     * @return bool
     */
    friend bool operator<(const UBInt &a, const UBInt &b);

    /**
     * @brief Tell if a is less than or equal to b.
     * @param a
     * @param b
     * @return bool
     */
    friend bool operator<=(const UBInt &a, const UBInt &b);

    /**
     * @brief Read a UBInt object via a std istream.
     * @return std::istream &
     */
    friend std::istream &operator>>(std::istream &, UBInt &);

    /**
     * @brief Print via a std ostream.
     * @return std::ostream &
     */
    friend std::ostream &operator<<(std::ostream &, const UBInt &);

    /**
     * @brief Transform to type u64. The result is not guarenteed if the input
     * is not a valid 64-bit integer.
     * @return u64
     */
    friend u64 to_u64(const UBInt &);

    /**
     * @brief Helper function: divide the input by 2.
     */
    friend void divide_by_2(UBInt &);

    /**
     * @brief Tell if the input is zero or not.
     * @return bool
     */
    friend bool is_zero(const UBInt &);

    /**
     * @brief Return the decimal length of this big integer.
     * @return int
     */
    friend int length(const UBInt &);

    /**
     * @brief Get the specified digit by index.
     * @param i Index of the needed digit.
     * @return int
     */
    int operator[](const int i) const;

private:
    /// A string storing all the digits as 8-bit integers (i.e. '\0' for digit
    /// 0, etc.).
    std::string digits_;
};

/**
 * @class CRTComposer
 * @brief Class for doing CRT composing on a specified group of moduli.
 *
 */
class CRTComposer {
public:
    /**
     * @brief Construct a new CRTComposer object from a set of moduli. Here the
     * CRT basis is produced.
     * @param moduli The CRT moduli.
     */
    CRTComposer(std::vector<u64> moduli);

    /**
     * @brief Compose a UBInt which is congruent to the remainders modulo
     * respective modulus.
     * @param remainders The remainders to compose.
     * @return UBInt
     */
    UBInt compose(std::vector<u64> remainders);

private:
    /**
     * @brief Helper function for computing the x^(-1) modulo a prime modulus.
     * @param x The input to inverse.
     * @param modulus The modulus which is assumed to be prime and of length <
     * 64bit.
     * @return UBInt
     */
    UBInt inv_mod_prime(const UBInt &x, const u64 modulus);

    /// The whole modulus which is the product of CRT moduli.
    UBInt whole_modulus;

    /// This vector stores the inverses of moduli products modulo corresponding
    /// modulus, multiplied by the product itself. In particular for the i-th
    /// modulus p_i, basis_inv_[i] is ((p_1*p_2*...*p_{i-1}*p_{i+1}*...)^(-1)
    /// mod p_i) * (p_1*p_2*...*p_{i-1}*p_{i+1}*...).
    std::vector<UBInt> basis_;

    /// The number of moduli in CRT basis.
    std::size_t basis_size_;
};

} // namespace hehub
