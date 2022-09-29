#include "catch2/catch.hpp"
#include "common/bigintpoly.h"
#include "common/mod_arith.h"
#include "common/permutation.h"
#include "common/sampling.h"
#include "primitives/ckks/ckks.h"
#include <iomanip>
#include <iostream>

using namespace hehub;

namespace hehub {
void fft_negacyclic_inplace(cc_double *coeffs, size_t log_poly_len,
                            bool inverse = false);
} // namespace hehub

TEST_CASE("fft", "[.]") {
    size_t logn = 8;
    size_t n = 1 << logn;
    const double eps = std::pow(2.0, -40);

    std::vector<cc_double> coeffs(n);
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, 1);
    for (auto &c : coeffs) {
        c = std::complex(distribution(generator), distribution(generator));
    }

    SECTION("FFT forward") {
        auto coeffs_copy(coeffs);

        auto substitute = [&](cc_double x_val) {
            cc_double sum = 0;
            cc_double x_pow = 1;
            for (size_t i = 0; i < n; i++) {
                sum += coeffs_copy[i] * x_pow;
                x_pow *= x_val;
            }
            return sum;
        };

        fft_negacyclic_inplace(coeffs.data(), logn);

        for (size_t i = 0; i < n; i = i * 3 + 1) {
            REQUIRE(abs(coeffs[i] - substitute(std::polar(
                                        1.0, (i * 2 + 1) * M_PI / n))) < eps);
        }
    }
    SECTION("FFT round trip") {
        auto coeffs_copy(coeffs);

        fft_negacyclic_inplace(coeffs.data(), logn);
        fft_negacyclic_inplace(coeffs.data(), logn, true);

        bool all_close = true;
        for (size_t i = 0; i < n; i = i * 3 + 1) {
            all_close = all_close && (abs(coeffs[i] - coeffs_copy[i]) < eps);
        }
        REQUIRE(all_close);
    }
}

TEST_CASE("ckks encoding") {
    size_t poly_len = 32;
    PolyDimensions pt_dim{
        poly_len,
        2,
        {144115188075593729, 144115188068319233}}; // ~114 bits in total

    auto data_size = poly_len / 2;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, 1);

    auto eps = std::pow(2, -35);
    auto check_all_close = [=](const auto &vec1, const auto &vec2) {
        for (size_t i = 0; i < data_size; i++) {
            if (abs(vec1.at(i) - vec2.at(i)) > eps) {
                return false;
            }
        }
        return true;
    };

    SECTION("small real") {
        double scaling_factor = std::pow(2.0, 40);
        std::vector<double> data(data_size);
        for (auto &d : data) {
            d = distribution(generator);
        }

        CkksPt pt = ckks::simd_encode(data, scaling_factor, pt_dim);
        auto data_recovered = ckks::simd_decode(pt);

        REQUIRE(data_recovered.size() == data.size());
        REQUIRE(check_all_close(data, data_recovered));
    }
    SECTION("big real") {
        double scaling_factor = std::pow(2.0, 80);
        std::vector<double> data(data_size);
        for (auto &d : data) {
            d = distribution(generator);
        }

        CkksPt pt = ckks::simd_encode(data, scaling_factor, pt_dim);
        auto data_recovered = ckks::simd_decode(pt);

        REQUIRE(data_recovered.size() == data.size());
        REQUIRE(check_all_close(data, data_recovered));
    }
    SECTION("small complex") {
        double scaling_factor = std::pow(2.0, 40);
        std::vector<cc_double> data(data_size);
        for (auto &d : data) {
            d = {distribution(generator), distribution(generator)};
        }

        CkksPt pt = ckks::simd_encode(data, scaling_factor, pt_dim);
        auto data_recovered = ckks::simd_decode<cc_double>(pt);

        REQUIRE(data_recovered.size() == data.size());
        REQUIRE(check_all_close(data, data_recovered));
    }
    SECTION("big complex") {
        double scaling_factor = std::pow(2.0, 80);
        std::vector<cc_double> data(data_size);
        for (auto &d : data) {
            d = {distribution(generator), distribution(generator)};
        }

        CkksPt pt = ckks::simd_encode(data, scaling_factor, pt_dim);
        auto data_recovered = ckks::simd_decode<cc_double>(pt);

        REQUIRE(data_recovered.size() == data.size());
        REQUIRE(check_all_close(data, data_recovered));
    }
}

TEST_CASE("ckks rescaling") {
    size_t poly_len = 8;
    std::vector<u64> moduli{17179672577, 17179410433,
                            17176854529}; // ~34 bits each
    PolyDimensions ct_dim{poly_len, 3, moduli};

    CkksCt ct;
    ct.scaling_factor = std::pow(2.0, 80);
    for (auto &c : ct) {
        c = get_rand_uniform_poly(ct_dim, PolyRepForm::coeff);
    }

    // Compose the coefficients
    std::array composed{UBigIntPoly(ct[0]), UBigIntPoly(ct[1])};

    // Rescale
    for (auto &c : ct) {
        ntt_negacyclic_inplace_lazy(c);
    }
    ckks::rescale_inplace(ct);
    for (auto &c : ct) {
        intt_negacyclic_inplace_lazy(c);
        strict_reduce(c);
    }

    // Checks
    auto dropped_mod = *moduli.crbegin();
    auto half_dropped_mod = dropped_mod / 2;
    REQUIRE(ct[0].component_count() == 2);
    REQUIRE(ct[1].component_count() == 2);
    double eps = std::pow(2.0, -60);
    REQUIRE(abs(ct.scaling_factor - std::pow(2.0, 80) / dropped_mod) < eps);

    std::array composed_new{UBigIntPoly(ct[0]), UBigIntPoly(ct[1])};

    for (auto half : {0, 1}) {
        for (size_t i = 0; i < poly_len; i++) {
            REQUIRE((composed[half][i] + half_dropped_mod) / dropped_mod ==
                    composed_new[half][i]);
        }
    }
}
