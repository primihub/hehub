#include "bigint.h"
#include "rns.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>

namespace hehub {

void UBInt::normalize() {
    while (limbs_.size() > 1 && limbs_.back() == 0)
        limbs_.pop_back();
    if (limbs_.empty())
        limbs_.push_back(0);
}

int UBInt::compare(const UBInt &other) const {
    if (limbs_.size() != other.limbs_.size())
        return limbs_.size() < other.limbs_.size() ? -1 : 1;
    for (int i = (int)limbs_.size() - 1; i >= 0; i--) {
        if (limbs_[i] != other.limbs_[i])
            return limbs_[i] < other.limbs_[i] ? -1 : 1;
    }
    return 0;
}

UBInt::UBInt(u64 nr) { limbs_.push_back(nr); }

UBInt::UBInt(const std::string &str) : UBInt(0ULL) {
    for (char c : str) {
        if (!std::isdigit(c))
            throw std::invalid_argument("str containing non-digit.");
        mul_u64(10);
        limbs_[0] += (c - '0');
        // propagate carry if last addition overflowed
        if (limbs_[0] < (u64)(c - '0')) {
            // very rare: only when limbs_[0] overflowed, which means
            // limbs_[0] was 0xFFFFFFFFFFFFFFFF and adding (c - '0') overflowed
            // since mul_u64 can make limbs_[0] up to 0xFFFFFFFFFFFFFFFF - 9,
            // this is extremely unlikely but handle for correctness
            for (size_t i = 1;; i++) {
                if (i >= limbs_.size())
                    limbs_.push_back(0);
                limbs_[i]++;
                if (limbs_[i] != 0)
                    break;
            }
        }
    }
    normalize();
}



UBInt UBInt::from_double(const double d) {
    if (d < 0 || std::isnan(d) || std::isinf(d))
        throw std::invalid_argument("Negative or invalid input.");

    u64 bits;
    std::memcpy(&bits, &d, sizeof(bits));

    int64_t biased_exp = (bits >> 52) & 0x7FF;
    int64_t exp = biased_exp - 1023;
    u64 mant = (bits & 0x000FFFFFFFFFFFFFULL) | 0x0010000000000000ULL;

    if (exp < 0)
        return UBInt(0);

    UBInt result;
    result.limbs_[0] = mant;
    size_t shift = (size_t)exp - 52;
    size_t limb_shift = shift / 64;
    size_t bit_shift = shift % 64;

    if (limb_shift > 0) {
        result.limbs_.resize(result.limbs_.size() + limb_shift, 0);
        std::rotate(result.limbs_.rbegin(),
                    result.limbs_.rbegin() + limb_shift,
                    result.limbs_.rend());
        std::fill_n(result.limbs_.begin(), limb_shift, 0);
    }

    if (bit_shift > 0) {
        u64 carry = 0;
        for (auto &limb : result.limbs_) {
            u64 new_val = (limb << bit_shift) | carry;
            carry = limb >> (64 - bit_shift);
            limb = new_val;
        }
        if (carry)
            result.limbs_.push_back(carry);
    }

    result.normalize();
    return result;
}

UBInt &UBInt::operator++() {
    for (size_t i = 0;; i++) {
        if (i >= limbs_.size())
            limbs_.push_back(0);
        limbs_[i]++;
        if (limbs_[i] != 0)
            break;
    }
    normalize();
    return *this;
}

UBInt UBInt::operator++(int _dummy_) {
    UBInt aux = *this;
    ++(*this);
    return aux;
}

UBInt &UBInt::operator--() {
    if (limbs_.size() == 1 && limbs_[0] == 0)
        throw std::runtime_error("UNDERFLOW");
    for (size_t i = 0;; i++) {
        if (limbs_[i] > 0) {
            limbs_[i]--;
            break;
        }
        limbs_[i] = 0xFFFFFFFFFFFFFFFFULL;
    }
    normalize();
    return *this;
}

UBInt UBInt::operator--(int _dummy_) {
    UBInt aux = *this;
    --(*this);
    return aux;
}

UBInt &operator+=(UBInt &a, const UBInt &b) {
    size_t max_len = std::max(a.limbs_.size(), b.limbs_.size());
    a.limbs_.resize(max_len + 1, 0);

    u64 carry = 0;
    for (size_t i = 0; i < max_len; i++) {
        u64 bi = i < b.limbs_.size() ? b.limbs_[i] : 0;
        u64 s;
        bool c1 = __builtin_add_overflow(a.limbs_[i], bi, &s);
        bool c2 = __builtin_add_overflow(s, carry, &s);
        a.limbs_[i] = s;
        carry = c1 + c2;
    }
    a.limbs_[max_len] = carry;
    a.normalize();
    return a;
}

UBInt operator+(const UBInt &a, const UBInt &b) {
    UBInt temp = a;
    temp += b;
    return temp;
}

UBInt &operator-=(UBInt &a, const UBInt &b) {
    if (a < b)
        throw std::runtime_error("UNDERFLOW");

    u64 borrow = 0;
    for (size_t i = 0; i < a.limbs_.size(); i++) {
        u64 bi = i < b.limbs_.size() ? b.limbs_[i] : 0;
        u64 d;
        bool b1 = __builtin_sub_overflow(a.limbs_[i], bi, &d);
        bool b2 = __builtin_sub_overflow(d, borrow, &d);
        a.limbs_[i] = d;
        borrow = b1 + b2;
    }
    a.normalize();
    return a;
}

UBInt operator-(const UBInt &a, const UBInt &b) {
    UBInt temp = a;
    temp -= b;
    return temp;
}

void UBInt::mul_u64(const u64 x) {
    if (x == 0) {
        limbs_.assign(1, 0);
        return;
    }
    u64 carry = 0;
    for (auto &limb : limbs_) {
        u128 product = (u128)limb * x + carry;
        limb = (u64)product;
        carry = product >> 64;
    }
    if (carry)
        limbs_.push_back(carry);
}

UBInt &operator*=(UBInt &a, const UBInt &b) {
    if (is_zero(a) || is_zero(b)) {
        a = UBInt();
        return a;
    }
    if (b.limbs_.size() == 1) {
        a.mul_u64(b.limbs_[0]);
        return a;
    }
    if (a.limbs_.size() == 1) {
        // a is single limb, b is multi — swap for efficiency
        u64 a0 = a.limbs_[0];
        a = b;
        a.mul_u64(a0);
        return a;
    }

    size_t n = a.limbs_.size(), m = b.limbs_.size();
    std::vector<u64> result(n + m, 0);

    for (size_t i = 0; i < n; i++) {
        u64 carry = 0;
        for (size_t j = 0; j < m; j++) {
            u128 product =
                (u128)a.limbs_[i] * b.limbs_[j] + result[i + j] + carry;
            result[i + j] = (u64)product;
            carry = product >> 64;
        }
        result[i + m] = carry;
    }

    a.limbs_ = std::move(result);
    a.normalize();
    return a;
}

UBInt operator*(const UBInt &a, const UBInt &b) {
    UBInt temp = a;
    temp *= b;
    return temp;
}

u64 UBInt::div_u64(const u64 divisor) {
    u64 remainder = 0;
    for (int i = (int)limbs_.size() - 1; i >= 0; i--) {
        u128 temp = ((u128)remainder << 64) | limbs_[i];
        limbs_[i] = (u64)(temp / divisor);
        remainder = (u64)(temp % divisor);
    }
    normalize();
    return remainder;
}

UBInt &operator/=(UBInt &a, const UBInt &b) {
    if (is_zero(b))
        throw std::invalid_argument("Arithmetic Error: Division By 0");
    int cmp = a.compare(b);
    if (cmp < 0) {
        a = UBInt();
        return a;
    }
    if (cmp == 0) {
        a = UBInt(1);
        return a;
    }
    if (b.limbs_.size() == 1) {
        a.div_u64(b.limbs_[0]);
        return a;
    }

    // Binary long division for multi-limb divisor
    size_t total_bits =
        a.limbs_.size() * 64 - (a.limbs_.back() ? __builtin_clzll(a.limbs_.back()) : 0);

    UBInt quotient(0);
    quotient.limbs_.resize(total_bits / 64 + 1, 0);

    UBInt remainder(0);

    for (int i = (int)total_bits - 1; i >= 0; i--) {
        // remainder <<= 1
        u64 carry = 0;
        for (auto &limb : remainder.limbs_) {
            u64 new_val = (limb << 1) | carry;
            carry = limb >> 63;
            limb = new_val;
        }
        if (carry)
            remainder.limbs_.push_back(1);

        // set bit i of a as LSB of remainder
        size_t limb_idx = (size_t)i / 64;
        size_t bit_idx = (size_t)i % 64;
        remainder.limbs_[0] |= (a.limbs_[limb_idx] >> bit_idx) & 1;

        if (remainder.compare(b) >= 0) {
            remainder -= b;
            size_t qi = (size_t)i / 64;
            size_t qb = (size_t)i % 64;
            quotient.limbs_[qi] |= (1ULL << qb);
        }
    }

    quotient.normalize();
    a = std::move(quotient);
    return a;
}

UBInt operator/(const UBInt &a, const UBInt &b) {
    UBInt temp = a;
    temp /= b;
    return temp;
}

UBInt &operator%=(UBInt &a, const UBInt &b) {
    a = a - (a / b) * b;
    return a;
}

UBInt operator%(const UBInt &a, const UBInt &b) {
    UBInt temp = a;
    temp %= b;
    return temp;
}

UBInt sqrt(const UBInt &a) {
    UBInt left(1), right(a), v(1), mid, prod;
    divide_by_2(right);
    while (left <= right) {
        mid += left;
        mid += right;
        divide_by_2(mid);
        prod = (mid * mid);
        if (prod <= a) {
            v = mid;
            ++mid;
            left = mid;
        } else {
            --mid;
            right = mid;
        }
        mid = UBInt();
    }
    return v;
}

bool operator==(const UBInt &a, const UBInt &b) {
    return a.compare(b) == 0;
}

bool operator!=(const UBInt &a, const UBInt &b) { return !(a == b); }

bool operator<(const UBInt &a, const UBInt &b) { return a.compare(b) < 0; }

bool operator>(const UBInt &a, const UBInt &b) { return a.compare(b) > 0; }

bool operator>=(const UBInt &a, const UBInt &b) { return a.compare(b) >= 0; }

bool operator<=(const UBInt &a, const UBInt &b) { return a.compare(b) <= 0; }

std::istream &operator>>(std::istream &in, UBInt &a) {
    std::string s;
    in >> s;
    a = UBInt(0);
    for (char c : s) {
        if (!std::isdigit(c))
            throw std::runtime_error("INVALID NUMBER");
        a.mul_u64(10);
        a.limbs_[0] += (c - '0');
        // propagate carry if overflow
        if (a.limbs_[0] < (u64)(c - '0')) {
            for (size_t i = 1;; i++) {
                if (i >= a.limbs_.size())
                    a.limbs_.push_back(0);
                a.limbs_[i]++;
                if (a.limbs_[i] != 0)
                    break;
            }
        }
    }
    a.normalize();
    return in;
}

std::ostream &operator<<(std::ostream &out, const UBInt &a) {
    if (is_zero(a)) {
        out << '0';
        return out;
    }
    // Convert to decimal by repeated division by 10
    std::string digits;
    UBInt temp(a);
    while (!is_zero(temp)) {
        u64 rem = temp.div_u64(10);
        digits.push_back('0' + (char)rem);
    }
    for (auto it = digits.rbegin(); it != digits.rend(); ++it)
        out << *it;
    return out;
}

u64 to_u64(const UBInt &a) {
    if (a.limbs_.empty())
        return 0;
    return a.limbs_[0];
}

double to_double(const UBInt &a) {
    if (is_zero(a))
        return 0.0;

    size_t top_idx = a.limbs_.size() - 1;
    u64 top = a.limbs_[top_idx];
    if (top == 0) {
        // shouldn't happen after normalize, but just in case
        return 0.0;
    }
    int lz = __builtin_clzll(top);
    int64_t top_bit_pos = (int64_t)top_idx * 64 + 63 - lz;

    int bits_from_top = 63 - lz;
    u64 mant;
    if (bits_from_top >= 52) {
        mant = (top & ((1ULL << bits_from_top) - 1)) >> (bits_from_top - 52);
    } else {
        u64 top_bits = top & ((1ULL << bits_from_top) - 1);
        int needed = 52 - bits_from_top;
        if (top_idx > 0) {
            u64 next = a.limbs_[top_idx - 1];
            top_bits = (top_bits << needed) | (next >> (64 - needed));
        } else {
            top_bits <<= needed;
        }
        mant = top_bits;
    }

    mant &= 0x000FFFFFFFFFFFFFULL;
    int64_t exp = top_bit_pos;
    u64 double_bits = ((u64)(exp + 1023) << 52) | mant;

    double result;
    std::memcpy(&result, &double_bits, sizeof(result));
    return result;
}

void divide_by_2(UBInt &a) {
    u64 carry = 0;
    for (int i = (int)a.limbs_.size() - 1; i >= 0; i--) {
        u64 new_val = (a.limbs_[i] >> 1) | (carry << 63);
        carry = a.limbs_[i] & 1;
        a.limbs_[i] = new_val;
    }
    a.normalize();
}

bool is_zero(const UBInt &a) {
    return a.limbs_.size() == 1 && a.limbs_[0] == 0;
}

int length(const UBInt &a) { return (int)a.limbs_.size(); }

int UBInt::operator[](const int index) const {
    if (index < 0)
        throw std::invalid_argument("Negative index");
    std::ostringstream oss;
    oss << *this;
    std::string s = oss.str();
    if ((size_t)index >= s.length())
        throw std::invalid_argument("Index out of range.");
    return s[s.length() - 1 - index] - '0';
}

CRTComposer::CRTComposer(std::vector<u64> moduli) {
    if (moduli.empty()) {
        throw std::invalid_argument("Empty CRT moduli set.");
    }
    for (size_t i = 0; i < moduli.size(); i++) {
        for (size_t j = i + 1; j < moduli.size(); j++) {
            if (moduli[j] == moduli[i]) {
                throw std::invalid_argument("Invalid CRT moduli set.");
            }
        }
    }

    basis_size_ = moduli.size();
    whole_modulus = UBInt(u64(1));
    for (size_t i = 0; i < basis_size_; i++) {
        whole_modulus *= UBInt(moduli[i]);
    }

    basis_.resize(basis_size_);
    for (size_t i = 0; i < basis_size_; i++) {
        UBInt prod(u64(1));
        for (size_t j = 0; j < basis_size_; j++) {
            if (j == i) {
                continue;
            }
            prod *= UBInt(moduli[j]);
        }
        UBInt inv = inv_mod_prime(prod, moduli[i]);
        basis_[i] = inv * prod;
    }
}

UBInt CRTComposer::compose(std::vector<u64> remainders) {
    if (remainders.size() != basis_size_) {
        throw std::invalid_argument("Number of remainders doesn't match.");
    }
    UBInt result(u64(0));
    for (size_t i = 0; i < basis_size_; i++) {
        result += basis_[i] * UBInt(remainders[i]);
    }
    result %= whole_modulus;
    return result;
}

UBInt CRTComposer::inv_mod_prime(const UBInt &x, const u64 modulus) {
    const auto modulus_big_int = UBInt(modulus);
    const auto index = modulus - 2;
    auto x_power = UBInt(u64(1));
    auto mask = u64(1) << 63;
    while (mask > index) {
        mask >>= 1;
    }
    while (mask) {
        x_power *= x_power;
        x_power %= modulus_big_int;
        if (mask & index) {
            x_power *= x;
            x_power %= modulus_big_int;
        }
        mask >>= 1;
    }
    return x_power;
}

UBIntVec::UBIntVec(const RnsPolynomial &rns_poly) {
    const auto dimension(rns_poly.dimension());
    const auto component_count(rns_poly.component_count());
    CRTComposer crt_composer(rns_poly.modulus_vec());
    for (size_t i = 0; i < dimension; i++) {
        std::vector<u64> remainder_coeffs;
        for (size_t j = 0; j < component_count; j++) {
            remainder_coeffs.push_back(rns_poly[j][i]);
        }
        coeffs_.push_back(crt_composer.compose(remainder_coeffs));
    }
}

std::ostream &operator<<(std::ostream &out, const UBIntVec &big_int_poly) {
    for (size_t i = big_int_poly.coeffs_.size() - 1; i > 0; i--) {
        out << big_int_poly.coeffs_[i];
        out << "*X^" << i << " + ";
    }
    if (!big_int_poly.coeffs_.empty()) {
        out << big_int_poly.coeffs_[0];
    }
    return out;
}

} // namespace hehub
