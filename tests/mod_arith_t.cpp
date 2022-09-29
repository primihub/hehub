#include "catch2/catch.hpp"
#include "common/mod_arith.h"

using namespace hehub;

TEST_CASE("batched barrett") {
    const u64 modulus =
        GENERATE(65537, 33333333, 777777777777777, 1234567890111111111);
    const size_t vec_len = 1000;

    u64 seed = 42;
    u64 vec[vec_len], vec_copy[vec_len];
    for (size_t i = 0; i < vec_len; i++) {
        seed = (seed ^ 893758435427369) * 65536 + 945738773644543;
        vec[i] = seed;
        vec_copy[i] = seed % modulus;
    }

    batched_barrett_lazy(modulus, vec_len, vec);
    auto check_range = [=](bool status, u64 item) {
        return status && (item < 2 * modulus);
    };
    REQUIRE(std::accumulate(vec, vec + vec_len, true, check_range));
    u64 diff[vec_len];
    for (size_t i = 0; i < vec_len; i++) {
        diff[i] = vec[i] - vec_copy[i];
    }
    auto check_equal = [=](bool status, u64 item) {
        return status && (item == 0 || item == modulus);
    };
    REQUIRE(std::accumulate(diff, diff + vec_len, true, check_equal));
}

TEST_CASE("batched mul mod") {
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

    SECTION("batched_mul_mod_hybrid") {
        batched_mul_mod_hybrid(modulus, vec_len, f, g, h);
        for (size_t i = 1; i < vec_len; i *= 3) {
            REQUIRE(h[i] == (u128)f[i] * g[i] % modulus);
        }
    }
    SECTION("batched_mul_mod_barrett") {
        batched_mul_mod_barrett(modulus, vec_len, f, g, h);
        for (size_t i = 1; i < vec_len; i *= 3) {
            REQUIRE(h[i] == (u128)f[i] * g[i] % modulus);
        }
    }
}

TEST_CASE("montgomery") {
    const u64 modulus = 38589379749438777;
    const size_t vec_len = 8;

    u128 seed = 42;
    std::vector<u128> f(vec_len);
    std::vector<u64> f_reduced(vec_len);
    for (size_t i = 0; i < vec_len; i++) {
        seed = seed * 3405898573857435 + 4385837453598385;
        f[i] = seed % ((u128)modulus << 64);
    }

    batched_montgomery_128_lazy(modulus, vec_len, f.data(), f_reduced.data());
    for (size_t i = 0; i < vec_len; i++) {
        REQUIRE(f_reduced[i] < 2 * modulus);
        REQUIRE(((u128)1 << 64) * f_reduced[i] % modulus == f[i] % modulus);
    }
}
