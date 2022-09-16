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
