#include "catch2/catch.hpp"
#include "common/bigintpoly.h"
#include "common/mod_arith.h"
#include "common/permutation.h"
#include "common/sampling.h"
#include "primitives/ckks/ckks.h"
#include <type_traits>

using namespace hehub;

namespace hehub {
void fft_negacyclic_natural_inout(cc_double *coeffs, size_t log_poly_len,
                                  bool inverse = false);
} // namespace hehub

#define REQUIRE_ALL_CLOSE(vec1, vec2, eps)                                     \
    do {                                                                       \
        static_assert(std::is_same<decltype(vec1), decltype(vec2)>::value);    \
        REQUIRE((vec1).size() == (vec2).size());                               \
        for (size_t i = 0; i < (vec1).size(); i++) {                           \
            REQUIRE(std::abs((vec1)[i] - (vec2)[i]) < (eps));                  \
        }                                                                      \
    } while (false);

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

        fft_negacyclic_natural_inout(coeffs.data(), logn);

        for (size_t i = 0; i < n; i = i * 3 + 1) {
            REQUIRE(abs(coeffs[i] - substitute(std::polar(
                                        1.0, (i * 2 + 1) * M_PI / n))) < eps);
        }
    }
    SECTION("FFT round trip") {
        auto coeffs_copy(coeffs);

        fft_negacyclic_natural_inout(coeffs.data(), logn);
        fft_negacyclic_natural_inout(coeffs.data(), logn, true);

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

    SECTION("small real") {
        double scaling_factor = std::pow(2.0, 40);
        std::vector<double> data(data_size);
        for (auto &d : data) {
            d = distribution(generator);
        }

        CkksPt pt = ckks::simd_encode(data, scaling_factor, pt_dim);
        auto data_recovered = ckks::simd_decode(pt);

        REQUIRE(data_recovered.size() == data.size());
        REQUIRE_ALL_CLOSE(data, data_recovered, pow(2.0, -35));
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
        REQUIRE_ALL_CLOSE(data, data_recovered, pow(2.0, -35));
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
        REQUIRE_ALL_CLOSE(data, data_recovered, pow(2.0, -35));
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
        REQUIRE_ALL_CLOSE(data, data_recovered, pow(2.0, -35));
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

TEST_CASE("ckks encryption") {
    std::vector<u64> ct_moduli{1099510054913}; // 40-bit
    size_t poly_len = 8;
    PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};
    RlweSk sk(ct_poly_dim);

    // random ckks plain data
    auto data_count = poly_len / 2;
    std::vector<double> plain_data(data_count);
    std::default_random_engine generator;
    std::normal_distribution<double> data_dist(0, 1);
    for (auto &d : plain_data) {
        d = data_dist(generator);
    }

    // encode
    int scaling_bits = 30;
    double scaling_factor = std::pow(2.0, scaling_bits);
    auto pt = ckks::simd_encode(plain_data, scaling_factor, ct_poly_dim);

    // encrypt
    auto ct = ckks::encrypt(pt, sk);

    // decrypt & decode
    auto data_recovered = ckks::simd_decode(ckks::decrypt(ct, sk));

    // check
    REQUIRE(data_recovered.size() == poly_len / 2);
    double eps = std::pow(2.0, 5 - scaling_bits); // noise's 6σ = 19.2 < 2^5
    REQUIRE_ALL_CLOSE(plain_data, data_recovered, eps);
}

TEST_CASE("ckks arith") {
    std::vector<u64> ct_moduli{1099510054913, 1073479681, 1072496641};
    size_t poly_len = 8;
    PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};
    RlweSk sk(ct_poly_dim);

    // random ckks plain data
    auto data_count = poly_len / 2;
    std::vector<double> plain_data1(data_count);
    std::vector<double> plain_data2(data_count);
    std::vector<double> plain_data3(data_count);
    std::default_random_engine generator;
    std::normal_distribution<double> data_dist(0, 1);
    for (auto &d : plain_data1) {
        d = data_dist(generator);
    }
    for (auto &d : plain_data2) {
        d = data_dist(generator);
    }
    for (auto &d : plain_data3) {
        d = data_dist(generator);
    }

    // encode
    int scaling_bits = 30;
    double scaling_factor = std::pow(2.0, scaling_bits);
    auto pt1 = ckks::simd_encode(plain_data1, scaling_factor, ct_poly_dim);
    auto pt2 = ckks::simd_encode(plain_data2, scaling_factor, ct_poly_dim);
    auto pt3 = ckks::simd_encode(plain_data3, scaling_factor, ct_poly_dim);

    SECTION("addition") {
        auto data_sum(plain_data1);
        for (size_t i = 0; i < data_count; i++) {
            data_sum[i] += plain_data2[i] + plain_data3[i];
        }

        // encrypt
        auto ct1 = ckks::encrypt(pt1, sk);
        auto ct2 = ckks::encrypt(pt2, sk);

        // add
        auto ct_sum = ckks::add(ct1, ct2);
        ct_sum = ckks::add_plain(ct_sum, pt3);

        // decrypt & decode
        auto sum_recovered = ckks::simd_decode(ckks::decrypt(ct_sum, sk));

        // check
        double eps = std::pow(2.0, 5 + 1 - scaling_bits);
        REQUIRE_ALL_CLOSE(data_sum, sum_recovered, eps);
    }
    SECTION("subtraction") {
        auto data_diff(plain_data1);
        for (size_t i = 0; i < data_count; i++) {
            data_diff[i] -= plain_data2[i] + plain_data3[i];
        }

        // encrypt
        auto ct1 = ckks::encrypt(pt1, sk);
        auto ct2 = ckks::encrypt(pt2, sk);

        // sub
        auto ct_diff = ckks::sub(ct1, ct2);
        ct_diff = ckks::sub_plain(ct_diff, pt3);

        // decrypt & decode
        auto diff_recovered = ckks::simd_decode(ckks::decrypt(ct_diff, sk));

        // check
        double eps = std::pow(2.0, 5 + 1 - scaling_bits);
        REQUIRE_ALL_CLOSE(data_diff, diff_recovered, eps);
    }
    SECTION("multiplication with plaintext") {
        auto data_prod(plain_data1);
        for (size_t i = 0; i < data_count; i++) {
            data_prod[i] *= plain_data2[i];
        }

        // encrypt
        auto ct1 = ckks::encrypt(pt1, sk);

        // mult
        auto ct_prod = ckks::mult_plain(ct1, pt2);

        SECTION("without rescaling") {
            // decrypt & decode
            auto prod_recovered = ckks::simd_decode(ckks::decrypt(ct_prod, sk));

            // check
            double eps =
                pow(2, 3 + 5 - scaling_bits); // abs of data < 6σ
                                              // with σ = data's std dev
            REQUIRE_ALL_CLOSE(data_prod, prod_recovered, eps);
        }
        SECTION("with rescaling") {
            // rescaling
            ckks::rescale_inplace(ct_prod);

            // decrypt & decode
            auto prod_recovered = ckks::simd_decode(ckks::decrypt(ct_prod, sk));

            // check
            double eps =
                pow(2, 3 + 5 - scaling_bits); // abs of data < 6σ
                                              // with σ = data's std dev
            REQUIRE_ALL_CLOSE(data_prod, prod_recovered, eps);
        }
    }
    SECTION("multiplication with ciphertext") {
        auto data_prod(plain_data1);
        for (size_t i = 0; i < data_count; i++) {
            data_prod[i] *= plain_data2[i];
        }

        // encrypt
        auto ct1 = ckks::encrypt(pt1, sk);
        auto ct2 = ckks::encrypt(pt2, sk);

        // gen relin key
        auto additional_mod = 1099507695617;
        auto relin_key = get_relin_key(sk, additional_mod);

        // mult
        auto ct_prod_quadratic = ckks::mult_low_level(ct1, ct2);
        auto ct_prod = ckks::relinearize(ct_prod_quadratic, relin_key);

        SECTION("without rescaling") {
            // decrypt & decode
            auto prod_recovered = ckks::simd_decode(ckks::decrypt(ct_prod, sk));

            // check
            double eps =
                pow(2, 3 + 5 + 1 - scaling_bits); // abs of data < 6σ
                                                  // with σ = data's std dev
            REQUIRE_ALL_CLOSE(data_prod, prod_recovered, eps);
        }
        SECTION("with rescaling") {
            // rescaling
            ckks::rescale_inplace(ct_prod);
            REQUIRE(ct_prod[0].component_count() == 2);
            REQUIRE(ct_prod[1].component_count() == 2);

            // decrypt & decode
            auto prod_recovered = ckks::simd_decode(ckks::decrypt(ct_prod, sk));

            // check
            double eps =
                pow(2, 3 + 5 + 1 - scaling_bits); // abs of data < 6σ
                                                  // with σ = data's std dev
            REQUIRE_ALL_CLOSE(data_prod, prod_recovered, eps);
        }
    }
}

TEST_CASE("ckks key switch") {
    SECTION("general key switching") {
        std::vector<u64> ct_moduli{1099510054913, 1073479681, 1072496641};
        size_t poly_len = 8;
        PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};
        u64 additional_mod = 1099507695617;

        auto data_count = poly_len / 2;
        std::vector<double> plain_data(data_count);
        std::default_random_engine generator;
        std::normal_distribution<double> data_dist(0, 1);
        for (auto &d : plain_data) {
            d = data_dist(generator);
        }

        RlweSk sk1(ct_poly_dim);
        RlweSk sk2(ct_poly_dim);
        auto ksk = RlweKsk(sk1, sk2, additional_mod);

        auto pt = ckks::simd_encode(plain_data, pow(2.0, 30), ct_poly_dim);
        auto ct = ckks::encrypt(pt, sk1);
        CkksCt ct_new = ext_prod_montgomery(ct[1], ksk);
        ckks::rescale_inplace(ct_new);
        ct_new.scaling_factor = ct.scaling_factor;
        ct_new[0] += ct[0];

        auto pt_recovered = ckks::decrypt(ct_new, sk2);
        auto data_recovered = ckks::simd_decode(pt_recovered);
        double eps = pow(2.0, -25); // empirical, needs analysis
        REQUIRE_ALL_CLOSE(plain_data, data_recovered, eps);
    }
    SECTION("conjugation") {
        std::vector<u64> ct_moduli{1099510054913, 1073479681, 1072496641};
        size_t poly_len = 8;
        PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};
        u64 additional_mod = 1099507695617;

        auto data_count = poly_len / 2;
        std::vector<cc_double> plain_data(data_count);
        std::vector<cc_double> data_conj;
        std::default_random_engine generator;
        std::normal_distribution<double> data_dist(0, 1);
        for (auto &d : plain_data) {
            d = {data_dist(generator), data_dist(generator)};
            data_conj.push_back(std::conj(d));
        }

        RlweSk sk(ct_poly_dim);
        auto conj_key = get_conj_key(sk, additional_mod);

        auto pt = ckks::simd_encode(plain_data, pow(2.0, 30), ct_poly_dim);
        auto ct = ckks::encrypt(pt, sk);
        auto ct_conj = ckks::conjugate(ct, conj_key);

        auto pt_recovered = ckks::decrypt(ct_conj, sk);
        auto data_recovered = ckks::simd_decode<cc_double>(pt_recovered);
        double eps = pow(2.0, -24); // empirical, needs analysis
        REQUIRE_ALL_CLOSE(data_conj, data_recovered, eps);
    }
    SECTION("rotation") {
        std::vector<u64> ct_moduli{1099510054913, 1073479681, 1072496641};
        size_t poly_len = 8;
        PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};
        u64 additional_mod = 1099507695617;

        auto data_count = poly_len / 2;
        std::vector<cc_double> plain_data(data_count);
        std::vector<cc_double> data_rotated(data_count);
        std::default_random_engine generator;
        std::normal_distribution<double> data_dist(0, 1);
        for (auto &d : plain_data) {
            d = {data_dist(generator), data_dist(generator)};
        }
        size_t step = GENERATE(1, 2, 3);
        for (size_t i = 0; i < data_count; i++) {
            data_rotated[(i + step) % data_count] = plain_data[i];
        }

        RlweSk sk(ct_poly_dim);
        auto rot_key_for_the_step = get_rot_key(sk, additional_mod, step);

        auto pt = ckks::simd_encode(plain_data, pow(2.0, 30), ct_poly_dim);
        auto ct = ckks::encrypt(pt, sk);
        auto ct_conj = ckks::rotate(ct, rot_key_for_the_step, step);

        auto pt_recovered = ckks::decrypt(ct_conj, sk);
        auto data_recovered = ckks::simd_decode<cc_double>(pt_recovered);
        double eps = pow(2.0, -23); // empirical, needs analysis
        REQUIRE_ALL_CLOSE(data_rotated, data_recovered, eps);
    }
}
