#include "catch2/catch.hpp"
#include "common/mod_arith.h"
#include "primitives/bgv/bgv.h"

using namespace hehub;

TEST_CASE("bgv encoding") {
    u64 p = 65537;
    size_t n = 1024;
    std::vector<u64> data(n);
    for (size_t i = 0; i < n; i++) {
        data[i] = (i * 888 + 123) % p;
    }

    auto pt = bgv::simd_encode(data, p);
    auto data_decoded = bgv::simd_decode(pt);
    REQUIRE(data_decoded == data);
}

TEST_CASE("bgv encryption") {
    SECTION("simple decrypt") {
        std::vector<u64> ct_moduli{35184358850561ULL, 36028796997599233ULL};
        u64 pt_modulus = 65537;
        size_t poly_len = 4096;
        PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};
        PolyDimensions pt_poly_dim{poly_len, 1, std::vector{pt_modulus}};

        RlweSk sk(ct_poly_dim);
        RlwePt pt(pt_poly_dim);
        // plaintext data
        u64 seed = 42;
        auto &pt_poly = pt[0];
        for (auto &coeff : pt_poly) {
            seed = (seed ^ 943598839) * 1024 + 348398497;
            coeff = seed % pt_modulus;
        }

        // encrypt & decrypt
        RlweCt ct = bgv::encrypt(pt, sk);
        RlwePt pt_recovered = bgv::decrypt(ct, sk, pt_modulus);

        // check
        REQUIRE(pt == pt_recovered);
    }
    SECTION("CRT decrypt") {
        std::vector<u64> ct_moduli{35184358850561ULL,
                                   36028796997599233ULL};
        size_t poly_len = 8;
        u64 pt_modulus = 576460752272228353ULL;

        PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};
        PolyDimensions pt_poly_dim{poly_len, 1, std::vector{pt_modulus}};

        RlweSk sk(ct_poly_dim);
        RlwePt pt(pt_poly_dim);
        // plaintext data
        u64 seed = 42;
        auto &pt_poly = pt[0];
        for (auto &coeff : pt_poly) {
            seed = (seed ^ 943598839) * 1023 + 348398497;
            coeff = seed % pt_modulus;
        }

        // encrypt & decrypt
        RlweCt ct = bgv::encrypt(pt, sk);
        RlwePt pt_recovered = bgv::decrypt(ct, sk, pt_modulus);

        // check
        REQUIRE(pt == pt_recovered);
    }
    SECTION("coprime condition") {
        std::vector<u64> ct_moduli{65537};
        u64 pt_modulus = 65537;
        size_t poly_len = 1024;

        PolyDimensions ct_poly_dim{poly_len, ct_moduli.size(), ct_moduli};
        PolyDimensions pt_poly_dim{poly_len, 1, std::vector{pt_modulus}};

        RlweSk sk(ct_poly_dim);
        RlwePt pt(pt_poly_dim);

        REQUIRE_THROWS(bgv::encrypt(pt, sk));
    }
}

TEST_CASE("bgv mod switch") {
    std::vector<u64> ct_moduli{131530753, 130809857};
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

    // original decryption result
    RlwePt pt = bgv::decrypt(ct, sk, pt_modulus);

    // mod switch and new decryption result
    bgv::mod_switch_inplace(ct, pt_modulus);
    RlwePt pt_new = bgv::decrypt(ct, sk, pt_modulus);

    // check
    bool all_equal = true;
    auto dropped_q = *ct_moduli.rbegin();
    for (size_t i = 0; i < poly_len; i++) {
        if ((u128)pt_new[0][i] * dropped_q % pt_modulus != pt[0][i]) {
            all_equal = false;
        }
    }
    REQUIRE(all_equal);
}
