#include "mod_arith.h"
#include <map>

namespace hehub {

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

inline u64 get_neg_q_inv(const u64 modulus) {
    auto tmp = std::get<1>(xgcd128((i128)1 << 64, ((i128)1 << 64) - modulus));
    return (u64)(tmp + ((i128)1 << 64));
}

inline u64 get_mont(const u64 modulus) {
    const u64 mont_but_one = (u64)(-1LL) % modulus;
    return mont_but_one + 1;
}

inline u64 get_mont_harvay(const u64 modulus) {
    const u64 mont_but_one = (u64)(-1LL) % modulus;
    return (((u128)(mont_but_one + 1)) << 64) / modulus;
}

void vector_mul_mod_hybrid_lazy(const u64 modulus, const size_t vec_len,
                                const u64 in_vec1[], const u64 in_vec2[],
                                u64 out_vec[]) {
    static std::map<u64, MulModLUT> lut_cache;

    auto it = lut_cache.find(modulus);
    if (it == lut_cache.end()) {
        lut_cache.insert(std::make_pair(
            modulus, std::make_tuple(get_neg_q_inv(modulus), get_mont(modulus),
                                     get_mont_harvay(modulus))));
        it = lut_cache.find(modulus);
    }
    const auto &lut = it->second;

    const u64 mshift = 64;
    const u64 neg_qinv = std::get<0>(lut);
    const u64 mont = std::get<1>(lut);
    const u64 mont_harvay = std::get<2>(lut);

    for (size_t i = 0; i < vec_len; i++) {
        u128 a = (u128)in_vec1[i] * in_vec2[i];
        u128 u = a * neg_qinv;
        u &= (((u128)1 << mshift) - 1);
        u *= modulus;

        u64 out_temp = (a + u) >> mshift;
        u64 out_temp2 = (u128)out_temp * mont_harvay >> 64;
        out_vec[i] = (u128)out_temp * mont - (u128)out_temp2 * modulus;
    }
}

void vector_mul_mod_barrett_lazy(const u64 modulus, const size_t vec_len,
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

void strict_reduce(RnsPolynomial &rns_poly) {
    auto mod_ptr = rns_poly.moduli_vec().begin();

    for (auto &component_poly : rns_poly) {
        auto curr_mod = *(mod_ptr++);
        for (auto &coeff : component_poly) {
            // here "coeff" can also be NTT value
            coeff -= (coeff >= curr_mod) ? curr_mod : 0;
        }
    }
}

} // namespace hehub
