#include "ckks.h"
#include "common/bigint.h"
#include "common/bigintpoly.h"
#include "common/mod_arith.h"
#include "common/rns_transform.h"
#include <numeric>

using namespace std;

namespace hehub {

void fft_negacyclic_inplace(cc_double *coeffs, size_t log_poly_len,
                            bool inverse = false) {
    auto poly_len = 1ULL << log_poly_len;
    vector<cc_double> coeffs_copy(poly_len);
    if (!inverse) {
        for (size_t i = 0; i < poly_len; i++) {
            coeffs_copy[i] = coeffs[i] * polar(1.0, i * M_PI / poly_len);
        }
    } else {
        for (size_t i = 0; i < poly_len; i++) {
            coeffs_copy[i] = coeffs[i];
        }
    }

    cc_double zeta = polar(1.0, 2 * M_PI / poly_len);
    if (inverse) {
        zeta = complex{zeta.real(), -zeta.imag()};
    }
    vector<cc_double> zpow{1.0};
    for (size_t i = 1; i < poly_len; i++) {
        zpow.push_back(zpow[i - 1] * zeta);
    }

    auto subs_zpow = [&](auto index) {
        cc_double sum = 0;
        for (size_t i = 0; i < poly_len; i++) {
            sum += coeffs_copy[i] * zpow[index * i % poly_len];
        }
        return sum;
    };

    vector<cc_double> values(poly_len);
    for (size_t i = 0; i < poly_len; i++) {
        values[i] = subs_zpow(i);
    }

    for (size_t i = 0; i < poly_len; i++) {
        coeffs[i] = values[i];
    }
    if (inverse) {
        for (size_t i = 0; i < poly_len; i++) {
            coeffs[i] *= polar(1.0, i * M_PI / poly_len * -1.0);
            coeffs[i] /= poly_len;
        }
    }
}

CkksPt ckks::simd_encode_cc(const vector<cc_double> &data,
                            const double scaling_factor,
                            const PolyDimensions &pt_poly_dim) {
    if (scaling_factor <= 0) {
        throw invalid_argument("Scaling factor should be positive.");
    }
    auto slot_count = pt_poly_dim.poly_len / 2;
    auto data_size = data.size();
    if (data_size > slot_count) {
        throw invalid_argument("Cannot encode " + to_string(data_size) +
                               " data into " + to_string(slot_count) +
                               " slots.");
    }

    vector<cc_double> interpolated(slot_count * 2, 0.0);
    for (size_t i = 0; i < data.size(); i++) {
        interpolated[i] = data[i];
        interpolated[slot_count * 2 - 1 - i] = conj(data[i]);
    }

    // Interpolate the data into element in C[X]/(X^n+1).
    fft_negacyclic_inplace(interpolated.data(),
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

    return pt;
}

vector<cc_double> ckks::simd_decode_cc(const CkksPt &pt,
                                       const double scaling_factor,
                                       size_t data_size) {
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
    vector<cc_double> data(pt.poly_len());

    auto pt_reduced(pt);
    strict_reduce(pt_reduced);
    auto poly_len = pt.poly_len();
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

    if (small_coeff) {
        auto first_mod = pt_reduced.modulus_at(0);
        auto half_first_mod = first_mod / 2;
        for (size_t i = 0; i < poly_len; i++) {
            if (pt_reduced[0][i] < half_first_mod) {
                data[i] = (double)pt_reduced[0][i];
            } else {
                data[i] = -(double)(first_mod - pt_reduced[0][i]);
            }
        }
    } else {
        UBigIntPoly pt_poly_big_int(pt_reduced);
        auto whole_modulus =
            accumulate(pt.modulus_vec().begin(), pt.modulus_vec().end(),
                       UBInt(1), [](auto acc, auto x) { return acc * x; });
        auto half_whole_mod = whole_modulus / 2;
        for (size_t i = 0; i < poly_len; i++) {
            stringstream ss;
            double abs_real;
            if (pt_poly_big_int[i] < half_whole_mod) {
                ss << pt_poly_big_int[i];
                ss >> abs_real;
                data[i] = abs_real;
            } else {
                ss << whole_modulus - pt_poly_big_int[i];
                ss >> abs_real;
                data[i] = -abs_real;
            }
        }
    }

    // Recover the data by scaling back and FFT
    for (auto &d : data) {
        d /= scaling_factor;
    }
    // currently data.size() == poly_len
    fft_negacyclic_inplace(data.data(), size_t(log2(data.size())));

    // Resize while removing conjugate part
    data.resize(data_size);
    return data;
}

template <>
vector<cc_double> ckks::simd_decode(const CkksPt &pt,
                                    const double scaling_factor,
                                    size_t data_size) {
    return simd_decode_cc(pt, scaling_factor, data_size);
}

template <>
vector<double> ckks::simd_decode(const CkksPt &pt, const double scaling_factor,
                                 size_t data_size) {
    auto data_cc = simd_decode_cc(pt, scaling_factor, data_size);
    vector<double> data;
    for (const auto d : data_cc) {
        data.push_back(d.real());
    }
    return data;
}

} // namespace hehub
