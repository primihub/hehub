#include "bigint.h"
#include "rns.h"
#include <sstream>

namespace hehub {

UBInt::UBInt(u64 nr) {
    do {
        digits_.push_back(nr % 10);
        nr /= 10;
    } while (nr);
}

UBInt::UBInt(const std::string &str) {
    digits_ = "";
    int n = str.size();
    for (int i = n - 1; i >= 0; i--) {
        if (!isdigit(str[i])) {
            throw std::invalid_argument("str containing non-digit.");
        }
        digits_.push_back(str[i] - '0');
    }
}

UBInt::UBInt(const char *str) {
    digits_ = "";
    for (int i = strlen(str) - 1; i >= 0; i--) {
        if (!isdigit(str[i])) {
            throw std::invalid_argument("str containing non-digit.");
        }
        digits_.push_back(str[i] - '0');
    }
}

UBInt UBInt::from_double(const double d) {
    if (d < 0) {
        throw std::invalid_argument("Negative input.");
    }

    std::stringstream ss;
    ss << std::fixed << d;
    std::string str;
    ss >> str;
    str.erase(str.begin() + str.find('.'), str.end());
    return UBInt(str);
}

UBInt &UBInt::operator++() {
    int i, n = digits_.size();
    for (i = 0; i < n && digits_[i] == 9; i++)
        digits_[i] = 0;
    if (i == n)
        digits_.push_back(1);
    else
        digits_[i]++;
    return *this;
}

UBInt UBInt::operator++(int _dummy_) {
    UBInt aux;
    aux = *this;
    ++(*this);
    return aux;
}

UBInt &UBInt::operator--() {
    if (digits_[0] == 0 && digits_.size() == 1)
        throw("UNDERFLOW");
    int i, n = digits_.size();
    for (i = 0; digits_[i] == 0 && i < n; i++)
        digits_[i] = 9;
    digits_[i]--;
    if (n > 1 && digits_[n - 1] == 0)
        digits_.pop_back();
    return *this;
}

UBInt UBInt::operator--(int _dummy_) {
    UBInt aux;
    aux = *this;
    --(*this);
    return aux;
}

UBInt &operator+=(UBInt &a, const UBInt &b) {
    int t = 0, s, i;
    int n = length(a), m = length(b);
    if (m > n)
        a.digits_.append(m - n, 0);
    n = length(a);
    for (i = 0; i < n; i++) {
        if (i < m)
            s = (a.digits_[i] + b.digits_[i]) + t;
        else
            s = a.digits_[i] + t;
        t = s / 10;
        a.digits_[i] = (s % 10);
    }
    if (t)
        a.digits_.push_back(t);
    return a;
}

UBInt operator+(const UBInt &a, const UBInt &b) {
    UBInt temp;
    temp = a;
    temp += b;
    return temp;
}

UBInt &operator-=(UBInt &a, const UBInt &b) {
    if (a < b)
        throw("UNDERFLOW");
    int n = length(a), m = length(b);
    int i, t = 0, s;
    for (i = 0; i < n; i++) {
        if (i < m)
            s = a.digits_[i] - b.digits_[i] + t;
        else
            s = a.digits_[i] + t;
        if (s < 0)
            s += 10, t = -1;
        else
            t = 0;
        a.digits_[i] = s;
    }
    while (n > 1 && a.digits_[n - 1] == 0)
        a.digits_.pop_back(), n--;
    return a;
}

UBInt operator-(const UBInt &a, const UBInt &b) {
    UBInt temp;
    temp = a;
    temp -= b;
    return temp;
}

UBInt &operator*=(UBInt &a, const UBInt &b) {
    if (is_zero(a) || is_zero(b)) {
        a = UBInt();
        return a;
    }
    int n = a.digits_.size(), m = b.digits_.size();
    std::vector<int> v(n + m, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            v[i + j] += (a.digits_[i]) * (b.digits_[j]);
        }
    }
    n += m;
    a.digits_.resize(v.size());
    for (int s, i = 0, t = 0; i < n; i++) {
        s = t + v[i];
        v[i] = s % 10;
        t = s / 10;
        a.digits_[i] = v[i];
    }
    for (int i = n - 1; i >= 1 && !v[i]; i--)
        a.digits_.pop_back();
    return a;
}

UBInt operator*(const UBInt &a, const UBInt &b) {
    UBInt temp;
    temp = a;
    temp *= b;
    return temp;
}

UBInt &operator/=(UBInt &a, const UBInt &b) {
    if (is_zero(b)) {
        throw std::invalid_argument("Arithmetic Error: Division By 0");
    }
    if (a < b) {
        a = UBInt();
        return a;
    }
    if (a == b) {
        a = UBInt(1);
        return a;
    }
    int i, lgcat = 0, cc;
    int n = length(a), m = length(b);
    std::vector<int> cat(n, 0);
    UBInt t;
    for (i = n - 1; t * 10 + a.digits_[i] < b; i--) {
        t *= 10;
        t += a.digits_[i];
    }
    for (; i >= 0; i--) {
        t = t * 10 + a.digits_[i];
        for (cc = 9; cc * b > t; cc--)
            ;
        t -= cc * b;
        cat[lgcat++] = cc;
    }
    a.digits_.resize(cat.size());
    for (i = 0; i < lgcat; i++)
        a.digits_[i] = cat[lgcat - i - 1];
    a.digits_.resize(lgcat);
    return a;
}

UBInt operator/(const UBInt &a, const UBInt &b) {
    UBInt temp;
    temp = a;
    temp /= b;
    return temp;
}

UBInt &operator%=(UBInt &a, const UBInt &b) {
    a -= a / b * b;
    return a;
}

UBInt operator%(const UBInt &a, const UBInt &b) {
    UBInt temp;
    temp = a;
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
    return a.digits_ == b.digits_;
}

bool operator!=(const UBInt &a, const UBInt &b) { return !(a == b); }

bool operator<(const UBInt &a, const UBInt &b) {
    int n = length(a), m = length(b);
    if (n != m)
        return n < m;
    while (n--)
        if (a.digits_[n] != b.digits_[n])
            return a.digits_[n] < b.digits_[n];
    return false;
}

bool operator>(const UBInt &a, const UBInt &b) { return b < a; }

bool operator>=(const UBInt &a, const UBInt &b) { return !(a < b); }

bool operator<=(const UBInt &a, const UBInt &b) { return !(a > b); }

std::istream &operator>>(std::istream &in, UBInt &a) {
    std::string s;
    in >> s;
    int n = s.size();
    for (int i = n - 1; i >= 0; i--) {
        if (!isdigit(s[i])) {
            throw std::runtime_error("INVALID NUMBER");
        }
        a.digits_[n - i - 1] = s[i];
    }
    return in;
}

std::ostream &operator<<(std::ostream &out, const UBInt &a) {
    for (int i = a.digits_.size() - 1; i >= 0; i--)
        out << (short)a.digits_[i];
    return out;
}

u64 to_u64(const UBInt &a) {
    std::stringstream ss;
    ss.str("");
    ss << a;
    u64 result;
    ss >> result;
    return result;
}

double to_double(const UBInt &a) {
    std::stringstream ss;
    ss.str("");
    ss << a;
    double result;
    ss >> result;
    return result;
}

void divide_by_2(UBInt &a) {
    int add = 0;
    for (int i = a.digits_.size() - 1; i >= 0; i--) {
        int digit = (a.digits_[i] >> 1) + add;
        add = ((a.digits_[i] & 1) * 5);
        a.digits_[i] = digit;
    }
    while (a.digits_.size() > 1 && !a.digits_.back()) {
        a.digits_.pop_back();
    }
}

bool is_zero(const UBInt &a) {
    if (a.digits_.size() == 1 && a.digits_[0] == 0)
        return true;
    return false;
}

int length(const UBInt &a) { return a.digits_.size(); }

int UBInt::operator[](const int index) const {
    if (digits_.size() <= index || index < 0)
        throw("ERROR");
    return digits_[index];
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
    for (size_t i = big_int_poly.coeffs_.size() - 1; i >= 0; i--) {
        out << big_int_poly.coeffs_[i];
        if (i == 0) { break; }
        out << "*X^" << i << " + ";
    }
    return out;
}

} // namespace hehub
