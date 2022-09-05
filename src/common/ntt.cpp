#include "ntt.h"
#include "mod_arith.h"
#include <cmath>
#include <map>

namespace hehub {

u64 __pow_mod(u64 modulus, u64 base, size_t index) {
    u64 power = 1;
    size_t mask = 1;
    while (mask <= index) {
        mask <<= 1;
    }
    mask >>= 1;
    while (mask) {
        power = (u128)power * power % modulus;
        if (mask & index) {
            power = (u128)power * base % modulus;
        }
        mask >>= 1;
    }
    return power;
}

u64 __find_a_2nth_unity_root(u64 modulus, u64 n) {
    if ((modulus - 1) % (2 * n) != 0) {
        throw std::invalid_argument("2N doesn't divide (modulus - 1)");
    }
    u64 candidate = 2;
    while (true) {
        if (__pow_mod(modulus, candidate, (modulus - 1) / 2) == modulus - 1) {
            break;
        }
        candidate++;
    }
    auto root = __pow_mod(modulus, candidate, (modulus - 1) / (2 * n));
    return root;
}

u64 __reverse_bits(u64 x, const size_t bit_len) {
    u64 mh = (1 << (bit_len - 1));
    u64 ml = 1;
    for (; mh >= ml; mh >>= 1, ml <<= 1) {
        if ((mh & x) && (ml & x)) {
        } else if (!(mh & x) && !(ml & x)) {
        } else if (!(mh & x) && (ml & x)) {
            x |= mh;
            x ^= ml;
        } else {
            x |= ml;
            x ^= mh;
        }
    }
    return x;
}

struct NTTFactors {
    NTTFactors(u64 modulus, size_t log_poly_len, bool for_inverse = false) {
        const size_t log_modulus = (u64)(log2(modulus) + 0.5);
        if (log_modulus > 59) {
            throw std::invalid_argument(
                "NTT not supporting primes with bit size > 59 currently.");
        }
        size_t poly_len = 1 << log_poly_len;

        const u64 root_of_2nth = __find_a_2nth_unity_root(modulus, poly_len);
        if (!for_inverse) {
            seq.resize(poly_len);
            seq_harvey.resize(poly_len);
            for (size_t i = 0; i < poly_len; i++) {
                seq[i] = __pow_mod(modulus, root_of_2nth,
                                   __reverse_bits(i, log_poly_len));
                seq_harvey[i] = ((u128)seq[i] << 64) / modulus;
            }
        } else {
            seq.resize(poly_len * 2);
            seq_harvey.resize(poly_len * 2);
            auto root_of_2nth_inv =
                __pow_mod(modulus, root_of_2nth, 2 * poly_len - 1);
            for (size_t l = 0; l < log_poly_len; l++) {
                auto start = (1 << l) - 1;
                auto power_index_factor = 1 << (log_poly_len - l);
                for (size_t i = 0; i < (1 << l); i++) {
                    auto idx = start + i;
                    seq[idx] =
                        __pow_mod(modulus, root_of_2nth_inv,
                                  __reverse_bits(i, l) * power_index_factor);
                    seq_harvey[idx] = ((u128)seq[idx] << 64) / modulus;
                }
            }
            const u64 poly_len_inv = modulus - ((modulus - 1) >> log_poly_len);
            const u64 poly_len_inv_harvey =
                ((u128)poly_len_inv << 64) / modulus;
            for (size_t i = 0; i < poly_len; i++) {
                auto temp = __pow_mod(modulus, root_of_2nth_inv, i);
                auto idx = i + poly_len;
                seq[idx] = mul_mod_harvey_lazy(modulus, temp, poly_len_inv,
                                               poly_len_inv_harvey);
                seq[idx] -= (seq[idx] >= modulus) ? modulus : 0;
                seq_harvey[idx] = ((u128)seq[idx] << 64) / modulus;
            }

            shuffled_indices.resize(poly_len);
            for (size_t i = 0; i < poly_len; i++) {
                shuffled_indices[i] = __reverse_bits(i, log_poly_len);
            }
        }
    }

    NTTFactors(const NTTFactors &copying) = default;

    NTTFactors(NTTFactors &&moving) = default;

    ~NTTFactors() {}

    std::vector<u64> seq;

    std::vector<u64> seq_harvey;

    std::vector<size_t> shuffled_indices;
};

std::map<std::pair<u64, u64>, NTTFactors> &ntt_factors_cache() {
    static std::map<std::pair<u64, u64>, NTTFactors> global_ntt_factors_cache;
    return global_ntt_factors_cache;
}

std::map<std::pair<u64, u64>, NTTFactors> &intt_factors_cache() {
    static std::map<std::pair<u64, u64>, NTTFactors> global_intt_factors_cache;
    return global_intt_factors_cache;
}

inline const auto &
__find_or_create_ntt_factors(const u64 modulus, const size_t log_poly_len,
                             const bool for_inverse = false) {
    if (!for_inverse) {
        auto it =
            ntt_factors_cache().find(std::make_pair(modulus, log_poly_len));
        if (it == ntt_factors_cache().end()) {
            ntt_factors_cache().insert(
                std::make_pair(std::make_pair(modulus, log_poly_len),
                               NTTFactors(modulus, log_poly_len)));
            it =
                ntt_factors_cache().find(std::make_pair(modulus, log_poly_len));
        }
        return it->second;
    } else {
        auto it =
            intt_factors_cache().find(std::make_pair(modulus, log_poly_len));
        if (it == intt_factors_cache().end()) {
            intt_factors_cache().insert(
                std::make_pair(std::make_pair(modulus, log_poly_len),
                               NTTFactors(modulus, log_poly_len, true)));
            it = intt_factors_cache().find(
                std::make_pair(modulus, log_poly_len));
        }
        return it->second;
    }
}

void ntt_negacyclic_inplace_lazy(const size_t log_poly_len, const u64 modulus,
                                 u64 coeffs[]) {
    const size_t poly_len = 1ULL << log_poly_len;
    // generate or read from cache
    const auto &ntt_factors =
        __find_or_create_ntt_factors(modulus, log_poly_len);

    size_t level, start, data_step, h, l;
    size_t idx = 1;
    u64 temp, zeta, zeta_harvey;
    for (level = 1, data_step = poly_len; level <= log_poly_len;
         level++, data_step >>= 1) {
        auto gap = data_step / 2;
        for (start = 0; start < poly_len; start += data_step, idx++) {
            zeta = ntt_factors.seq[idx];
            zeta_harvey = ntt_factors.seq_harvey[idx];
            for (l = start; l < start + gap; l++) {
                h = l + gap;
                temp =
                    mul_mod_harvey_lazy(modulus, coeffs[h], zeta, zeta_harvey);
                coeffs[h] = coeffs[l] + 2 * modulus - temp;
                coeffs[l] = coeffs[l] + temp;
            }
        }
    }

    const u64 log_modulus = (u64)(log2(modulus) + 0.5);
    const u64 div_fix = (modulus >= (1ULL << log_modulus)) ? 1 : 0;
    for (size_t i = 0; i < poly_len; i++) {
        coeffs[i] -= ((coeffs[i] >> log_modulus) - div_fix) * modulus;
    }
}

void intt_negacyclic_inplace_lazy(const size_t log_poly_len, const u64 modulus,
                                  u64 values[]) {
    const size_t poly_len = 1ULL << log_poly_len;
    // generate or read from cache
    const auto &intt_factors =
        __find_or_create_ntt_factors(modulus, log_poly_len, true);

    u64 values_shuffled[poly_len];
    const auto &shuffled_indices = intt_factors.shuffled_indices;
    for (size_t i = 0; i < poly_len; i++) {
        values_shuffled[i] = values[shuffled_indices[i]];
    }

    size_t level, start, data_step, h, l;
    size_t idx = 0;
    u64 temp, zeta, zeta_harvey;
    for (level = 1, data_step = poly_len; level <= log_poly_len;
         level++, data_step >>= 1) {
        auto gap = data_step / 2;
        for (start = 0; start < poly_len; start += data_step, idx++) {
            zeta = intt_factors.seq[idx];
            zeta_harvey = intt_factors.seq_harvey[idx];
            for (l = start; l < start + gap; l++) {
                h = l + gap;
                temp = mul_mod_harvey_lazy(modulus, values_shuffled[h], zeta,
                                           zeta_harvey);
                values_shuffled[h] = values_shuffled[l] + 2 * modulus - temp;
                values_shuffled[l] = values_shuffled[l] + temp;
            }
        }
    }

    for (size_t i = 0; i < poly_len; i++) {
        values[i] = values_shuffled[shuffled_indices[i]];
    }

    const u64 log_modulus = (u64)(log2(modulus) + 0.5);
    const u64 div_fix = (modulus >= (1ULL << log_modulus)) ? 1 : 0;
    idx++;
    for (size_t i = 0; i < poly_len; i++, idx++) {
        values[i] -= ((values[i] >> log_modulus) - div_fix) * modulus;
        values[i] =
            mul_mod_harvey_lazy(modulus, values[i], intt_factors.seq[idx],
                                intt_factors.seq_harvey[idx]);
    }
}

void ntt_negacyclic_inplace_lazy(RnsPolynomial &rns_poly) {
    const auto component_count = rns_poly.components_.size();
    auto mod_ptr = rns_poly.moduli_.begin();
    const auto log_poly_len = rns_poly.log_poly_len_;

    for (auto &component_poly : rns_poly) {
        ntt_negacyclic_inplace_lazy(log_poly_len, *(mod_ptr++),
                                    component_poly.data());
    }

    rns_poly.rep_form = PolyRepForm::value;
}

void intt_negacyclic_inplace_lazy(RnsPolynomial &rns_poly) {
    const auto component_count = rns_poly.components_.size();
    auto mod_ptr = rns_poly.moduli_.begin();
    const auto log_poly_len = rns_poly.log_poly_len_;

    for (auto &component_poly : rns_poly) {
        intt_negacyclic_inplace_lazy(log_poly_len, *(mod_ptr++),
                                     component_poly.data());
    }

    rns_poly.rep_form = PolyRepForm::coeff;
}

void cache_ntt_factors_strict(const u64 log_poly_len,
                              const std::vector<u64> &moduli) {
    for (auto modulus : moduli) {
        __find_or_create_ntt_factors(modulus, log_poly_len);
        __find_or_create_ntt_factors(modulus, log_poly_len, true);
    }
}

} // namespace hehub
