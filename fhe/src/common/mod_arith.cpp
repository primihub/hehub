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

inline u64 get_mont(const u64 modulus) {
    const u64 mont_but_one = (u64)(-1LL) % modulus;
    return mont_but_one + 1;
}

inline u64 get_mont_harvey(const u64 modulus) {
    const u64 mont_but_one = (u64)(-1LL) % modulus;
    return (((u128)(mont_but_one + 1)) << 64) / modulus;
}

void batched_mul_mod_hybrid_lazy(const u64 modulus, const size_t vec_len,
                                 const u64 in_vec1[], const u64 in_vec2[],
                                 u64 out_vec[]) {
    static std::map<u64, MulModLUT> lut_cache;

    auto it = lut_cache.find(modulus);
    if (it == lut_cache.end()) {
        lut_cache.insert(std::make_pair(
            modulus,
            std::make_tuple(get_inv_neg_q_mod_2to64(modulus), get_mont(modulus),
                            get_mont_harvey(modulus))));
        it = lut_cache.find(modulus);
    }
    const auto [neg_qinv, mont, mont_harvey] = it->second;
    const u64 mshift = 64;

    for (size_t i = 0; i < vec_len; i++) {
        // The Montgomery part
        u128 a = (u128)in_vec1[i] * in_vec2[i];
        u128 u = a * neg_qinv;
        u &= (((u128)1 << mshift) - 1);
        u *= modulus;

        // The D. Harvey part
        u64 out_temp = (a + u) >> mshift;
        u64 out_temp2 = (u128)out_temp * mont_harvey >> 64;
        out_vec[i] = (u128)out_temp * mont - (u128)out_temp2 * modulus;
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
