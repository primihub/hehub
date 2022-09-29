#include "mod_arith.h"
#include <cmath>
#include <map>

namespace hehub {

using MulModLUT = std::tuple<u64, u64, u64>;

void batched_barrett_lazy(const u64 modulus, const size_t vec_len, u64 vec[]) {
    u64 c = (u64)(-1) / modulus;
    for (size_t i = 0; i < vec_len; i++) {
        u128 a = (u128)vec[i] * c;
        u64 approx_quotient = a >> 64;
        u128 approx_mod_multiple = modulus * approx_quotient;
        vec[i] -= approx_mod_multiple;
    }
}

std::tuple<i128, i128, i128> xgcd128(i128 a, i128 b) {
    int sign_a = (a < 0) ? (-1) : 1;
    a *= sign_a;
    int sign_b = (b < 0) ? (-1) : 1;
    b *= sign_b;

    i128 prev_x = 1;
    i128 x = 0;
    i128 prev_y = 0;
    i128 y = 1;
    i128 t = 0;
    while (b) {
        i128 q = a / b;

        t = prev_x - x * q;
        prev_x = x;
        x = t;

        t = prev_y - y * q;
        prev_y = y;
        y = t;

        t = b;
        b = a - b * q;
        a = t;
    }

    return std::make_tuple(sign_a * prev_x, sign_b * prev_y, a);
}

inline u64 get_inv_neg_q_mod_2to64(const u64 modulus) {
    auto tmp = std::get<1>(xgcd128((i128)1 << 64, ((i128)1 << 64) - modulus));
    return (u64)(tmp + ((i128)1 << 64));
}

inline u64 get_2to64_reduced(const u64 modulus) {
    const u64 _2to64_but_one = (u64)(-1LL) % modulus;
    return _2to64_but_one + 1;
}

inline u64 get_2to64_harvey(const u64 modulus) {
    const u64 _2to64_but_one_reduced = (u64)(-1LL) % modulus;
    return (((u128)(_2to64_but_one_reduced + 1)) << 64) / modulus;
}

void batched_mul_mod_hybrid_lazy(const u64 modulus, const size_t vec_len,
                                 const u64 in_vec1[], const u64 in_vec2[],
                                 u64 out_vec[]) {
    static std::map<u64, MulModLUT> lut_cache;

    auto it = lut_cache.find(modulus);
    if (it == lut_cache.end()) {
        lut_cache.insert(std::make_pair(
            modulus,
            std::make_tuple(get_inv_neg_q_mod_2to64(modulus), get_2to64_reduced(modulus),
                            get_2to64_harvey(modulus))));
        it = lut_cache.find(modulus);
    }
    const auto [neg_qinv, _2to64_reduced, _2to64_harvey] = it->second;
    const u64 mshift = 64;

    for (size_t i = 0; i < vec_len; i++) {
        // The Montgomery part
        u128 a = (u128)in_vec1[i] * in_vec2[i];
        u128 u = a * neg_qinv;
        u &= (((u128)1 << mshift) - 1);
        u *= modulus;

        // The D. Harvey part
        u64 out_temp = (a + u) >> mshift;
        u64 out_temp2 = (u128)out_temp * _2to64_harvey >> 64;
        out_vec[i] = (u128)out_temp * _2to64_reduced - (u128)out_temp2 * modulus;
    }
}

void batched_mul_mod_barrett_lazy(const u64 modulus, const size_t vec_len,
                                  const u64 in_vec1[], const u64 in_vec2[],
                                  u64 out_vec[]) {
    u128 c = (u128)(-1) / modulus;
    u64 ch = c >> 64;
    u64 cl = c;
    for (size_t i = 0; i < vec_len; i++) {
        u128 a = (u128)in_vec1[i] * in_vec2[i];
        u64 ah = a >> 64;
        u64 al = a;
        u64 ah_ch = ah * ch;
        u128 ah_cl = (u128)ah * cl;
        u128 al_ch = (u128)al * ch;
        u64 approx_quotient = ah_ch + ((ah_cl + al_ch) >> 64);
        u128 approx_mod_multiple = (u128)modulus * approx_quotient;
        out_vec[i] = a - approx_mod_multiple;
    }
}

void batched_montgomery_128_lazy(const u64 modulus, const size_t len,
                                 const u128 in[], u64 out[]) {
    static std::map<u64, u64> neg_q_inv_lut;
    u64 neg_q_inv;

    auto it = neg_q_inv_lut.find(modulus);
    if (it == neg_q_inv_lut.end()) {
        neg_q_inv = get_inv_neg_q_mod_2to64(modulus);
        neg_q_inv_lut.insert(std::make_pair(modulus, neg_q_inv));
    } else {
        neg_q_inv = it->second;
    }

    for (int i = 0; i < len; i++) {
        u128 a = in[i];
        u128 u = a * neg_q_inv;
        u = (u64)u;
        u *= modulus;
        a = (a + u) >> 64;
        out[i] = a;
    }
}

std::map<std::pair<u64, u64>, u64> modular_inverse_table;

u64 inverse_mod_prime(const u64 elem, const u64 prime) {
    auto it = modular_inverse_table.find(std::make_pair(elem, prime));
    if (it != modular_inverse_table.end()) {
        return it->second;
    } else {
        auto result = std::get<1>(xgcd128((i128)prime, (i128)elem));
        modular_inverse_table.insert(
            std::make_pair(std::make_pair(elem, prime), result));
        return result;
    }
}

} // namespace hehub
