#include "catch2/catch.hpp"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include "common/sampling.h"
#include "bgv/bgv.h"
#include <iostream>

using namespace hehub;

TEST_CASE("bgv encoding") {
    u64 p = 65537;
    size_t n = 1024;
    std::vector<u64> data(n);
    u64 seed = 1;
    for (auto &d : data) {
        d = ((seed++) * 888 + 123) % p;
    }

    // encode and decode
    auto pt = bgv::simd_encode(data, p);
    auto data_decoded = bgv::simd_decode(pt);

    // check
    REQUIRE(data_decoded == data);
}

TEST_CASE("bgv encryption") {
    std::vector<u64> ct_moduli{131530753, 130809857};
    size_t dimension = 128;
    RnsPolyParams ct_params{dimension, ct_moduli.size(), ct_moduli};

    RlweSk sk(ct_params);

    SECTION("simple decrypt") {
        u64 pt_modulus = 65537;
        RnsPolyParams pt_params{dimension, 1, std::vector{pt_modulus}};

        // random plaintext data
        BgvPt pt = get_rand_uniform_poly(pt_params);

        // encrypt & decrypt
        auto ct = bgv::encrypt(pt, sk);
        auto pt_recovered = bgv::decrypt(ct, sk);

        // check
        REQUIRE(pt == pt_recovered);
    }
    SECTION("CRT decrypt") {
        u64 pt_modulus = 35184358850561ULL;
        RnsPolyParams pt_params{dimension, 1, std::vector{pt_modulus}};

        // random plaintext data
        BgvPt pt = get_rand_uniform_poly(pt_params);

        // encrypt & decrypt
        auto ct = bgv::encrypt(pt, sk);
        auto pt_recovered = bgv::decrypt(ct, sk);

        // check
        REQUIRE(pt == pt_recovered);
    }
    SECTION("coprime condition") {
        // now plaintext modulus is not coprime with ciphertext modulus,
        // which is not allowed
        u64 pt_modulus = 131530753;
        RnsPolyParams pt_params{dimension, 1, std::vector{pt_modulus}};

        BgvPt pt(pt_params);

        REQUIRE_THROWS(bgv::encrypt(pt, sk));
    }
}

TEST_CASE("bgv arith") {
    std::vector<u64> ct_moduli{131530753, 130809857};
    size_t dimension = 8;
    RnsPolyParams ct_params{dimension, ct_moduli.size(), ct_moduli};

    RlweSk sk(ct_params);

    SECTION("addition") {
        u64 pt_modulus = 65537;
        RnsPolyParams pt_params{dimension, 1, std::vector{pt_modulus}};

        // random plaintext data
        BgvPt pt1 = get_rand_uniform_poly(pt_params);
        BgvPt pt2 = get_rand_uniform_poly(pt_params);
        BgvPt pt3 = get_rand_uniform_poly(pt_params);
        BgvPt pt_sum = pt1 + pt2 + pt3;
        reduce_strict(pt_sum);

        // encrypt
        auto ct1 = bgv::encrypt(pt1, sk);
        auto ct2 = bgv::encrypt(pt2, sk);

        // add
        auto ct_sum = bgv::add(ct1, ct2);
        ct_sum = bgv::add_plain(ct_sum, pt3);

        // decrypt
        auto sum_recovered = bgv::decrypt(ct_sum, sk);

        // check
        REQUIRE(pt_sum == sum_recovered);
    }
    SECTION("subtraction") {
        u64 pt_modulus = 65537;
        RnsPolyParams pt_params{dimension, 1, std::vector{pt_modulus}};

        // random plaintext data
        BgvPt pt1 = get_rand_uniform_poly(pt_params);
        BgvPt pt2 = get_rand_uniform_poly(pt_params);
        BgvPt pt3 = get_rand_uniform_poly(pt_params);
        BgvPt pt_diff = pt1 - pt2 - pt3;
        reduce_strict(pt_diff);

        // encrypt
        auto ct1 = bgv::encrypt(pt1, sk);
        auto ct2 = bgv::encrypt(pt2, sk);

        // sub
        auto ct_diff = bgv::sub(ct1, ct2);
        ct_diff = bgv::sub_plain(ct_diff, pt3);

        // decrypt
        auto diff_recovered = bgv::decrypt(ct_diff, sk);

        // check
        REQUIRE(pt_diff == diff_recovered);
    }
    SECTION("multiplication with plaintext") {
        u64 pt_modulus = 65537;
        RnsPolyParams pt_params{dimension, 1, std::vector{pt_modulus}};

        // random plaintext data
        BgvPt pt1 = get_rand_uniform_poly(pt_params);
        auto ct1 = bgv::encrypt(pt1, sk);
        BgvPt pt2 = get_rand_uniform_poly(pt_params);
        auto ct_prod = bgv::mult_plain(ct1, pt2);

        // mult on plaintexts
        ntt_negacyclic_inplace_lazy(pt1);
        ntt_negacyclic_inplace_lazy(pt2);
        BgvPt pt_prod = pt1 * pt2;
        intt_negacyclic_inplace_lazy(pt_prod);
        reduce_strict(pt_prod);

        // decrypt
        auto prod_recovered = bgv::decrypt(ct_prod, sk);

        // check
        REQUIRE(pt_prod == prod_recovered);
    }
    SECTION("encode & arith") {
        u64 pt_modulus = 65537;
        auto data_count = dimension;
        std::vector<u64> plain_data1(data_count);
        std::vector<u64> plain_data2(data_count);
        std::vector<u64> plain_data3(data_count);
        u64 seed = 1;
        for (auto &d : plain_data1) {
            d = ((seed++) * 888 + 123) % pt_modulus;
        }
        for (auto &d : plain_data2) {
            d = ((seed++) * 888 + 123) % pt_modulus;
        }
        for (auto &d : plain_data3) {
            d = ((seed++) * 888 + 123) % pt_modulus;
        }

        // encode
        auto pt1 = bgv::simd_encode(plain_data1, pt_modulus);
        auto pt2 = bgv::simd_encode(plain_data2, pt_modulus);
        auto pt3 = bgv::simd_encode(plain_data3, pt_modulus);

        // encrypt & arith
        auto ct1 = bgv::encrypt(pt1, sk);
        auto ct_prod = bgv::mult_plain(ct1, pt2);
        auto ct3 = bgv::encrypt(pt3, sk);
        auto ct_res = bgv::add(ct_prod, ct3);

        // decrypt
        auto res_decrypted = bgv::decrypt(ct_res, sk);
        auto res_data = bgv::simd_decode(res_decrypted);

        // check
        for (size_t i = 0; i < data_count; i++) {
            REQUIRE(res_data[i] ==
                    (plain_data1[i] * plain_data2[i] + plain_data3[i]) %
                        pt_modulus);
        }
    }
    // SECTION("multiplication with ciphertext") {
    //     u64 pt_modulus = 65537;
    //     auto data_count = dimension;
    //     std::vector<u64> plain_data1(data_count);
    //     std::vector<u64> plain_data2(data_count);
    //     u64 seed = 1;
    //     for (auto &d : plain_data1) {
    //         d = ((seed++) * 888 + 123) % pt_modulus;
    //     }
    //     for (auto &d : plain_data2) {
    //         d = ((seed++) * 888 + 123) % pt_modulus;
    //     }

    //     // encode
    //     auto pt1 = bgv::simd_encode(plain_data1, pt_modulus);
    //     auto pt2 = bgv::simd_encode(plain_data2, pt_modulus);

    //     // encrypt & arith
    //     auto ct1 = bgv::encrypt(pt1, sk);
    //     auto ct2 = bgv::encrypt(pt2, sk);
    //     auto ct_prod_quadratic = bgv::mult_low_level(ct1, ct2);
    //     auto relin_key = get_relin_key(sk, 131923969);
    //     auto ct_prod = bgv::relinearize(ct_prod_quadratic, relin_key);

    //     // decrypt
    //     auto prod_decrypted = bgv::decrypt(ct_prod, sk);
    //     auto prod_data = bgv::simd_decode(prod_decrypted);

    //     // check
    //     for (size_t i = 0; i < data_count; i++) {
    //         REQUIRE(prod_data[i] ==
    //                 plain_data1[i] * plain_data2[i] % pt_modulus);
    //     }
    // }
}

TEST_CASE("bgv mod switch") {
    std::vector<u64> ct_moduli{140737486520321, 140737485864961};
    u64 pt_modulus = 65537;
    size_t dimension = 8;
    RnsPolyParams ct_params{dimension, ct_moduli.size(), ct_moduli};
    RnsPolyParams pt_params{dimension, 1, std::vector{pt_modulus}};

    RlweSk sk(ct_params);
    RnsPolynomial fake_pt(ct_params);
    // random data, seen as real plaintext with large noise
    u64 seed = 42;
    for (size_t i = 0; i < dimension; i++) {
        seed = (seed * 1024 + 348398497) % 12345678901111111;
        for (size_t k = 0; k < ct_moduli.size(); k++) {
            fake_pt[k][i] = seed % ct_moduli[k];
        }
    }
    ntt_negacyclic_inplace_lazy(fake_pt);
    BgvCt ct = bgv::get_rlwe_sample_lift_noise(sk, pt_modulus);
    ct.plain_modulus = pt_modulus;
    ct[1] += fake_pt;

    // the "actual" plaintext
    auto pt = bgv::decrypt(ct, sk);

    // mod switch and new decryption result
    bgv::mod_switch_inplace(ct);
    auto pt_new = bgv::decrypt(ct, sk);

    REQUIRE(pt_new == pt);
}
