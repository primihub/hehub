#include "catch2/catch.hpp"
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

    std::vector<cc_double> coeffs(n);
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, 1);
    for (auto &c : coeffs) {
        c = std::complex(distribution(generator), distribution(generator));
    }
    auto coeffs_copy(coeffs);

    fft_negacyclic_inplace(coeffs.data(), logn);
    fft_negacyclic_inplace(coeffs.data(), logn, true);

    const double eps = std::pow(2.0, -45);
    bool all_close = true;
    for (size_t i = 0; i < n; i++) {
        if (abs(coeffs[i] - coeffs_copy[i]) > eps) {
            all_close = false;
            break;
        }
    }
    REQUIRE(all_close);
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

    SECTION("small complex") {
        double scaling_factor = std::pow(2.0, 40);
        std::vector<cc_double> data(data_size);
        for (auto &d : data) {
            d = {distribution(generator), distribution(generator)};
        }

        CkksPt pt = ckks::simd_encode(data, scaling_factor, pt_dim);
        auto data_recovered = ckks::simd_decode<cc_double>(pt, scaling_factor);

        REQUIRE(data_recovered.size() == data.size());
        REQUIRE(check_all_close(data, data_recovered));
    }
    SECTION("small real") {
        double scaling_factor = std::pow(2.0, 40);
        std::vector<double> data(data_size);
        for (auto &d : data) {
            d = distribution(generator);
        }

        CkksPt pt = ckks::simd_encode(data, scaling_factor, pt_dim);
        auto data_recovered = ckks::simd_decode(pt, scaling_factor);

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
        auto data_recovered = ckks::simd_decode<cc_double>(pt, scaling_factor);

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
        auto data_recovered = ckks::simd_decode(pt, scaling_factor);

        REQUIRE(data_recovered.size() == data.size());
        REQUIRE(check_all_close(data, data_recovered));
    }
}
