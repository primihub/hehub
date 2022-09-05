#include "catch2/catch.hpp"
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
        std::vector<u64> ct_moduli{35184358850561ULL, 36028796997599233ULL,
                                   576460752272228353ULL};
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
        std::vector<u64> ct_moduli{65537, 35184358850561ULL,
                                   36028796997599233ULL};
        size_t poly_len = 1024;
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
