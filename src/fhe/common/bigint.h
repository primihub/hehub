/**
 * @file ubigint.h
 * @brief Class of big integers and utility of CRT composing.
 *
 */

#pragma once

#include "type_defs.h"
#include "rns.h"
#include <cstring>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace hehub {

class UBInt {
public:
    UBInt(u64 nr = 0);

    UBInt(const std::string &str);

    UBInt(const char *str);

    static UBInt from_double(const double d);

    UBInt(const UBInt &other) = default;

    UBInt(UBInt &&other) = default;

    UBInt &operator=(const UBInt &other) = default;

    UBInt &operator=(UBInt &&other) = default;

    UBInt &operator++();

    UBInt operator++(int _dummy_);

    UBInt &operator--();

    UBInt operator--(int _dummy_);

    friend UBInt &operator+=(UBInt &a, const UBInt &b);

    friend UBInt operator+(const UBInt &a, const UBInt &b);

    friend UBInt &operator-=(UBInt &a, const UBInt &b);

    friend UBInt operator-(const UBInt &a, const UBInt &b);

    friend UBInt &operator*=(UBInt &a, const UBInt &b);

    friend UBInt operator*(const UBInt &a, const UBInt &b);

    friend UBInt &operator/=(UBInt &a, const UBInt &b);

    friend UBInt operator/(const UBInt &a, const UBInt &b);

    friend UBInt &operator%=(UBInt &a, const UBInt &b);

    friend UBInt operator%(const UBInt &a, const UBInt &b);

    friend UBInt sqrt(const UBInt &);

    friend bool operator==(const UBInt &a, const UBInt &b);

    friend bool operator!=(const UBInt &a, const UBInt &b);

    friend bool operator>(const UBInt &a, const UBInt &b);

    friend bool operator>=(const UBInt &a, const UBInt &b);

    friend bool operator<(const UBInt &a, const UBInt &b);

    friend bool operator<=(const UBInt &a, const UBInt &b);

    friend std::istream &operator>>(std::istream &, UBInt &);

    friend std::ostream &operator<<(std::ostream &, const UBInt &);

    friend u64 to_u64(const UBInt &);

    friend double to_double(const UBInt &);

    friend void divide_by_2(UBInt &);

    friend bool is_zero(const UBInt &);

    friend int length(const UBInt &);

    int operator[](const int i) const;

private:
    std::string digits_;
};

class CRTComposer {
public:
    CRTComposer(std::vector<u64> moduli);

    UBInt compose(std::vector<u64> remainders);

private:
    UBInt inv_mod_prime(const UBInt &x, const u64 modulus);

    UBInt whole_modulus;

    std::vector<UBInt> basis_;

    std::size_t basis_size_;
};

class UBIntVec {
public:

    UBIntVec(const RnsPolynomial &rns_poly);

    inline const size_t dimension() const { return coeffs_.size(); }

    inline UBInt &operator[](int i) { return coeffs_[i]; }

    inline const UBInt &operator[](int i) const { return coeffs_[i]; }

    friend std::ostream &operator<<(std::ostream &, const UBIntVec &);

private:
    std::vector<UBInt> coeffs_;
};

} // namespace hehub
