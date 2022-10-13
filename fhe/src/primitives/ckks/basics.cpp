#include "ckks.h"
#include "common/bigint.h"
#include "common/bigintpoly.h"
#include "common/mod_arith.h"
#include "common/permutation.h"
#include "common/rns_transform.h"
#include <numeric>

using namespace std;

namespace hehub {

/// @brief Inplace FFT with coefficients/point values input and output in
/// natural order.
void fft_negacyclic_natural_inout(cc_double *coeffs, size_t log_poly_len,
                                  bool inverse = false) {
    /**************************** Preparations ******************************/
    struct FFTFactors {
        FFTFactors(size_t log_poly_len, bool inverse = false) {
            auto poly_len = 1ULL << log_poly_len;
            cc_double zeta = polar(1.0, 2 * M_PI / poly_len);
            if (inverse) {
                zeta = conj(zeta);
            }
            for (size_t i = 0; i < poly_len; i++) {
                coeff_trans.push_back(
                    inverse ? polar(1.0 / poly_len, i * M_PI / poly_len * -1.0)
                            : polar(1.0, i * M_PI / poly_len));
            }

            size_t level, local_idx, gap;
            for (level = 1, gap = poly_len / 2; level <= log_poly_len;
                 level++, gap >>= 1) {
                for (local_idx = 0; local_idx < poly_len / gap / 2;
                     local_idx++) {
                    auto zeta_pow =
                        pow(zeta, (__bit_rev_naive_16(local_idx, level - 1)
                                   << (log_poly_len - level)));
                    butterfly.push_back(zeta_pow);
                }
            }
        }

        vector<cc_double> coeff_trans;

        vector<cc_double> butterfly;
    };

    static map<pair<size_t, bool>, FFTFactors> fft_factors_cache;

    auto __find_or_create_fft_factors = [&](size_t log_poly_len, bool inverse) {
        const auto args = make_pair(log_poly_len, inverse);
        auto it = fft_factors_cache.find(args);
        if (it == fft_factors_cache.end()) {
            fft_factors_cache.insert(
                std::make_pair(args, FFTFactors(log_poly_len, inverse)));
            it = fft_factors_cache.find(args);
        }
        return it->second;
    };

    const FFTFactors &fft_factors =
        __find_or_create_fft_factors(log_poly_len, inverse);
    /************************* End of preparations ***************************/

    auto poly_len = 1ULL << log_poly_len;
    vector<cc_double> coeffs_copy(poly_len);
    if (!inverse) {
        for (size_t i = 0; i < poly_len; i++) {
            coeffs_copy[i] = coeffs[i] * fft_factors.coeff_trans[i];
        }
    } else {
        for (size_t i = 0; i < poly_len; i++) {
            coeffs_copy[i] = coeffs[i];
        }
    }

    size_t level, start, gap, h, l, idx = 0;
    for (level = 1, gap = poly_len / 2; level <= log_poly_len;
         level++, gap >>= 1) {
        for (start = 0; start < poly_len; start += 2 * gap, idx++) {
            for (l = start; l < start + gap; l++) {
                h = l + gap;
                auto temp = coeffs_copy[h] * fft_factors.butterfly[idx];
                coeffs_copy[h] = coeffs_copy[l] - temp;
                coeffs_copy[l] = coeffs_copy[l] + temp;
            }
        }
    }

    for (size_t i = 0; i < poly_len; i++) {
        coeffs[i] = coeffs_copy[__bit_rev_naive_16(i, log_poly_len)];
    }
    if (inverse) {
        for (size_t i = 0; i < poly_len; i++) {
            coeffs[i] *= fft_factors.coeff_trans[i];
        }
    }
}

CkksPt ckks::simd_encode_cc(const vector<cc_double> &data,
                            const double scaling_factor,
                            const PolyDimensions &pt_poly_dim) {
    if (scaling_factor <= 0) {
        throw invalid_argument("Scaling factor should be positive.");
    }
    auto poly_len = pt_poly_dim.poly_len;
    size_t log_poly_len = round(log2(poly_len));
    auto slot_count = poly_len / 2;
    auto data_size = data.size();
    if (data_size > slot_count) {
        throw invalid_argument("Cannot encode " + to_string(data_size) +
                               " data into " + to_string(slot_count) +
                               " slots.");
    }

    vector<cc_double> interpolated(poly_len, 0.0);
    auto &root_indices = root_index_factors();
    auto mask = (1 << (log_poly_len + 1)) - 1; // for fast modulo 2*len
    for (size_t i = 0; i < data.size(); i++) {
        auto root_index = root_indices[i] & mask;
        auto position = (root_index - 1) / 2;
        interpolated[position] = data[i];
        interpolated[poly_len - 1 - position] = conj(data[i]);
    }

    // Interpolate the data into element in C[X]/(X^n+1).
    fft_negacyclic_natural_inout(interpolated.data(),
                                 size_t(log2(interpolated.size())),
                                 /*inverse=*/true);

    CkksPt pt(pt_poly_dim);

    // Decide if coefficients after scaling will be all smaller than 64-bit.
    bool small_coeff = true;
    double small_bound = pow(2.0, 64) / scaling_factor;
    for (const auto &d : interpolated) {
        if (abs(d.real()) > small_bound) {
            small_coeff = false;
            break;
        }
    }
    if (small_coeff) {
        // Transform the coefficients into u64.
        vector<u64> coeffs_u64_abs;
        vector<bool> coeffs_being_neg;
        for (auto &coeff_cc : interpolated) {
            double coeff_rr = coeff_cc.real();
            coeff_rr *= scaling_factor;
            coeffs_u64_abs.push_back(u64(abs(coeff_rr)));
            coeffs_being_neg.push_back(coeff_rr <= 0);
        }

        // Migrate the (abs values of) coefficients into RNS.
        for (size_t k = 0; k < pt_poly_dim.component_count; k++) {
            copy(coeffs_u64_abs.begin(), coeffs_u64_abs.end(), pt[k].data());
            batched_barrett_lazy(pt.modulus_at(k), pt.poly_len(), pt[k].data());
        }

        // Recover the signs of coefficients.
        for (size_t i = 0; i < pt.poly_len(); i++) {
            if (coeffs_being_neg[i]) {
                for (size_t k = 0; k < pt_poly_dim.component_count; k++) {
                    pt[k][i] =
                        (pt[k][i] == 0) ? 0 : (2 * pt.modulus_at(k) - pt[k][i]);
                }
            }
        }
    } else {
        // Transform the coefficients into big int.
        vector<UBInt> coeffs_bigint_abs;
        vector<bool> coeffs_being_neg;

        for (auto &coeff_cc : interpolated) {
            double coeff_rr = coeff_cc.real();
            coeff_rr *= scaling_factor;
            coeffs_being_neg.push_back(coeff_rr <= 0);
            coeffs_bigint_abs.push_back(UBInt::from_double(abs(coeff_rr)));
        }

        // Migrate the coefficients into RNS
        vector<UBInt> moduli_big;
        for (auto &mod : pt_poly_dim.moduli) {
            moduli_big.push_back(UBInt(mod));
        }
        for (size_t i = 0; i < pt.poly_len(); i++) {
            for (size_t k = 0; k < pt_poly_dim.component_count; k++) {
                pt[k][i] = to_u64(coeffs_bigint_abs[i] % moduli_big[k]);
            }

            // Recover the sign
            if (coeffs_being_neg[i]) {
                for (size_t k = 0; k < pt_poly_dim.component_count; k++) {
                    pt[k][i] = pt.modulus_at(k) - pt[k][i];
                }
            }
        }
    }

    pt.scaling_factor = scaling_factor;
    return pt;
}

vector<cc_double> ckks::simd_decode_cc(const CkksPt &pt, size_t data_size) {
    auto scaling_factor = pt.scaling_factor;
    if (scaling_factor <= 0) {
        throw invalid_argument("Scaling factor should be positive.");
    }
    auto slot_count = pt.poly_len() / 2;
    if (data_size == 0) {
        data_size = slot_count; // Actual default argument
    }
    if (data_size > slot_count) {
        throw invalid_argument("Cannot decode " + to_string(data_size) +
                               " items from " + to_string(slot_count) +
                               " slots.");
    }

    auto pt_reduced(pt);
    strict_reduce(pt_reduced);
    auto poly_len = pt.poly_len();
    size_t log_poly_len = round(log2(poly_len));
    auto components = pt.component_count();

    // Decide whether the coefficients when in composed form will be all smaller
    // than the first modulus.
    bool small_coeff = true;
    RnsPolynomial first_component(poly_len, 1, pt.modulus_vec());
    first_component[0] = pt_reduced[0];
    vector rest_moduli(pt.modulus_vec().begin() + 1, pt.modulus_vec().end());
    auto first_compo_under_rest_mod =
        rns_base_transform(first_component, rest_moduli);
    for (size_t k = 0; k < components - 1; k++) {
        if (first_compo_under_rest_mod[k] != pt_reduced[k + 1]) {
            small_coeff = false;
            break;
        }
    }

    vector<cc_double> interpolated(pt.poly_len());
    if (small_coeff) {
        auto first_mod = pt_reduced.modulus_at(0);
        auto half_first_mod = first_mod / 2;
        for (size_t i = 0; i < poly_len; i++) {
            if (pt_reduced[0][i] < half_first_mod) {
                interpolated[i] = (double)pt_reduced[0][i];
            } else {
                interpolated[i] = -(double)(first_mod - pt_reduced[0][i]);
            }
        }
    } else {
        UBigIntPoly pt_poly_big_int(pt_reduced);
        auto whole_modulus =
            accumulate(pt.modulus_vec().begin(), pt.modulus_vec().end(),
                       UBInt(1), [](auto acc, auto x) { return acc * x; });
        auto half_whole_mod = whole_modulus / 2;
        for (size_t i = 0; i < poly_len; i++) {
            double abs_real;
            if (pt_poly_big_int[i] < half_whole_mod) {
                abs_real = to_double(pt_poly_big_int[i]);
                interpolated[i] = abs_real;
            } else {
                abs_real = to_double(whole_modulus - pt_poly_big_int[i]);
                interpolated[i] = -abs_real;
            }
        }
    }

    // Recover the data by scaling back and FFT
    for (auto &i : interpolated) {
        i /= scaling_factor;
    }
    // interpolated.size() == poly_len
    fft_negacyclic_natural_inout(interpolated.data(), log_poly_len);

    // extract the original conjugation half
    vector<cc_double> data(slot_count);
    auto &root_indices = root_index_factors();
    auto mask = (1 << (log_poly_len + 1)) - 1; // for fast modulo 2*len
    for (size_t i = 0; i < slot_count; i++) {
        auto root_index = root_indices[i] & mask;
        auto position = (root_index - 1) / 2;
        data[i] = interpolated[position];
    }
    return data;
}

template <>
vector<cc_double> ckks::simd_decode(const CkksPt &pt, size_t data_size) {
    return simd_decode_cc(pt, data_size);
}

template <>
vector<double> ckks::simd_decode(const CkksPt &pt, size_t data_size) {
    auto data_cc = simd_decode_cc(pt, data_size);
    vector<double> data;
    for (const auto d : data_cc) {
        data.push_back(d.real());
    }
    return data;
}

} // namespace hehub
