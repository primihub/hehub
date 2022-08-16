#include "catch2/catch.hpp"
#include "common/mod_arith.h"

using namespace hehub;

TEST_CASE("vector mul mod") {
    const u64 modulus = 1234567890111111111;
    const size_t vec_len = 1000;

    u64 seed = 42;
    u64 f[vec_len], g[vec_len], h[vec_len];
    for (size_t i = 0; i < vec_len; i++) {
        seed = seed * 65968279837582827 ^ 3948528936546489545;
        f[i] = seed % modulus;
        seed = seed * 43534547657678213 ^ 7955436776934235466;
        g[i] = seed % modulus;
    }

    SECTION("vector_mul_mod_hybrid") {
        vector_mul_mod_hybrid(modulus, vec_len, f, g, h);
        for (size_t i = 1; i < vec_len; i *= 3) {
            REQUIRE(h[i] == (u128)f[i] * g[i] % modulus);
        }
    }
    SECTION("vector_mul_mod_barrett") {
        vector_mul_mod_barrett(modulus, vec_len, f, g, h);
        for (size_t i = 1; i < vec_len; i *= 3) {
            REQUIRE(h[i] == (u128)f[i] * g[i] % modulus);
        }
    }
}
