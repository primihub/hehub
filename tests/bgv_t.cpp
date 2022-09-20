#include "catch2/catch.hpp"
#include "common/mod_arith.h"
#include "common/sampling.h"
#include "primitives/bgv/bgv.h"

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
    size_t poly_len = 128;
    PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};

    RlweSk sk(ct_poly_dim);

    SECTION("simple decrypt") {
        u64 pt_modulus = 65537;
        PolyDimensions pt_poly_dim{poly_len, 1, std::vector{pt_modulus}};

        // random plaintext data
        RlwePt pt = get_rand_uniform_poly(pt_poly_dim);

        // encrypt & decrypt
        RlweCt ct = bgv::encrypt(pt, sk);
        RlwePt pt_recovered = bgv::decrypt(ct, sk, pt_modulus);

        // check
        REQUIRE(pt == pt_recovered);
    }
    SECTION("CRT decrypt") {
        u64 pt_modulus = 35184358850561ULL;
        PolyDimensions pt_poly_dim{poly_len, 1, std::vector{pt_modulus}};

        // random plaintext data
        RlwePt pt = get_rand_uniform_poly(pt_poly_dim);

        // encrypt & decrypt
        RlweCt ct = bgv::encrypt(pt, sk);
        RlwePt pt_recovered = bgv::decrypt(ct, sk, pt_modulus);

        // check
        REQUIRE(pt == pt_recovered);
    }
    SECTION("coprime condition") {
        // now plaintext modulus is not coprime with ciphertext modulus,
        // which is not allowed
        u64 pt_modulus = 131530753;
        PolyDimensions pt_poly_dim{poly_len, 1, std::vector{pt_modulus}};

        RlwePt pt(pt_poly_dim);

        REQUIRE_THROWS(bgv::encrypt(pt, sk));
    }
}

TEST_CASE("bgv arith") {
    std::vector<u64> ct_moduli{131530753, 130809857};
    size_t poly_len = 8;
    PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};

    RlweSk sk(ct_poly_dim);

    SECTION("addition") {
        u64 pt_modulus = 65537;
        PolyDimensions pt_poly_dim{poly_len, 1, std::vector{pt_modulus}};

        // random plaintext data
        RlwePt pt1 = get_rand_uniform_poly(pt_poly_dim);
        RlwePt pt2 = get_rand_uniform_poly(pt_poly_dim);
        RlwePt pt3 = get_rand_uniform_poly(pt_poly_dim);
        RlwePt pt_sum = pt1 + pt2 + pt3;
        strict_reduce(pt_sum);

        // encrypt
        auto ct1 = bgv::encrypt(pt1, sk);
        auto ct2 = bgv::encrypt(pt2, sk);

        // add
        auto ct_sum = bgv::add(ct1, ct2);
        ct_sum = bgv::add_plain(ct_sum, pt3);

        // decrypt
        RlwePt sum_recovered = bgv::decrypt(ct_sum, sk, pt_modulus);

        // check
        REQUIRE(pt_sum == sum_recovered);
    }
    SECTION("subtraction") {
        u64 pt_modulus = 65537;
        PolyDimensions pt_poly_dim{poly_len, 1, std::vector{pt_modulus}};

        // random plaintext data
        RlwePt pt1 = get_rand_uniform_poly(pt_poly_dim);
        RlwePt pt2 = get_rand_uniform_poly(pt_poly_dim);
        RlwePt pt3 = get_rand_uniform_poly(pt_poly_dim);
        RlwePt pt_diff = pt1 - pt2 - pt3;
        strict_reduce(pt_diff);

        // encrypt
        auto ct1 = bgv::encrypt(pt1, sk);
        auto ct2 = bgv::encrypt(pt2, sk);

        // sub
        auto ct_diff = bgv::sub(ct1, ct2);
        ct_diff = bgv::sub_plain(ct_diff, pt3);

        // decrypt
        RlwePt diff_recovered = bgv::decrypt(ct_diff, sk, pt_modulus);

        // check
        REQUIRE(pt_diff == diff_recovered);
    }
    SECTION("multiplication") {
        u64 pt_modulus = 65537;
        PolyDimensions pt_poly_dim{poly_len, 1, std::vector{pt_modulus}};

        // random plaintext data
        RlwePt pt1 = get_rand_uniform_poly(pt_poly_dim);
        auto ct1 = bgv::encrypt(pt1, sk);
        RlwePt pt2 = get_rand_uniform_poly(pt_poly_dim);
        auto ct_prod = bgv::mult_plain(ct1, pt2);

        // mult on plaintexts
        ntt_negacyclic_inplace_lazy(pt1);
        ntt_negacyclic_inplace_lazy(pt2);
        RlwePt pt_prod = pt1 * pt2;
        intt_negacyclic_inplace_lazy(pt_prod);
        strict_reduce(pt_prod);

        // decrypt
        RlwePt prod_recovered = bgv::decrypt(ct_prod, sk, pt_modulus);

        // check
        REQUIRE(pt_prod == prod_recovered);
    }
    SECTION("encode & arith") {
        u64 pt_modulus = 65537;
        auto data_count = poly_len;
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
        auto res_decrypted = bgv::decrypt(ct_res, sk, pt_modulus);
        auto res_data = bgv::simd_decode(res_decrypted);

        // check
        for (size_t i = 0; i < data_count; i = i * 2 + 1) {
            REQUIRE(res_data[i] ==
                    (plain_data1[i] * plain_data2[i] + plain_data3[i]) %
                        pt_modulus);
        }
    }
}

TEST_CASE("bgv mod switch") {
    std::vector<u64> ct_moduli{140737486520321, 140737485864961};
    u64 pt_modulus = 65537;
    size_t poly_len = 8;
    PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};
    PolyDimensions pt_poly_dim{poly_len, 1, std::vector{pt_modulus}};

    RlweSk sk(ct_poly_dim);
    RlwePt fake_pt(ct_poly_dim);
    // random data, seen as real plaintext with large noise
    u64 seed = 42;
    for (size_t i = 0; i < poly_len; i++) {
        seed = (seed * 1024 + 348398497) % 12345678901111111;
        for (size_t k = 0; k < ct_moduli.size(); k++) {
            fake_pt[k][i] = seed % ct_moduli[k];
        }
    }
    ntt_negacyclic_inplace_lazy(fake_pt);
    RlweCt ct = bgv::get_rlwe_sample_lift_noise(sk, ct_poly_dim, pt_modulus);
    ct[1] += fake_pt;

    // the "actual" plaintext
    RlwePt pt = bgv::decrypt(ct, sk, pt_modulus);

    // mod switch and new decryption result
    bgv::mod_switch_inplace(ct, pt_modulus);
    RlwePt pt_new = bgv::decrypt(ct, sk, pt_modulus);

    REQUIRE(pt_new == pt);
}
