#include "catch2/catch.hpp"
#include "fhe/common/mod_arith.h"

using namespace hehub;

TEST_CASE("BM reductions") {
    const u64 modulus = 1234567890111111111;
    const size_t vec_len = 8192;

    u64 seed = 42;
    u64 f[vec_len], g[vec_len], h[vec_len], k[vec_len];
    for (size_t i = 0; i < vec_len; i++) {
        seed = seed * 65968279837582827 ^ 3948528936546489545;
        f[i] = seed % modulus;
        seed = seed * 43534547657678213 ^ 7955436776934235466;
        g[i] = seed % modulus;
        seed = (seed ^ 39857467872338747) * 65536 + 394866313;
        h[i] = seed;
        seed = (seed ^ 39857467872338747) * 65536 + 394866313;
        k[i] = seed;
    }

    BENCHMARK("batched_barrett_lazy 8192") {
        return batched_barrett_lazy(modulus, vec_len, h);
    };

    BENCHMARK("batched_barrett 8192") {
        return batched_barrett(modulus, vec_len, k);
    };

    BENCHMARK("batched_mul_mod_hybrid_lazy 8192") {
        return batched_mul_mod_hybrid_lazy(modulus, vec_len, f, g, h);
    };

    BENCHMARK("batched_mul_mod_barrett_lazy 8192") {
        return batched_mul_mod_barrett_lazy(modulus, vec_len, f, g, h);
    };

    BENCHMARK("batched_mul_mod_hybrid 8192") {
        return batched_mul_mod_hybrid(modulus, vec_len, f, g, h);
    };

    BENCHMARK("batched_mul_mod_barrett 8192") {
        return batched_mul_mod_barrett(modulus, vec_len, f, g, h);
    };
}
