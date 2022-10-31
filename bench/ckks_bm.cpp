#include "catch2/catch.hpp"
#include "ckks/ckks.h"
#include "common/ntt.h"

using namespace hehub;
using namespace std;

TEST_CASE("benchmark ckks") {
    for (auto [LOGN, scaling_bits] :
         {pair{12, 36}, pair{13, 43}, pair{14, 48}, pair{15, 55}}) {
        u64 N = 1 << LOGN;
        auto params = ckks::create_params(N, scaling_bits);
        cache_ntt_factors_strict(LOGN, params.moduli);

        CkksSk sk(params);
        std::vector<cc_double> data(N / 2);
        BENCHMARK(string("CKKS encode+encrypt / N=") + to_string(N) +
                  string(" / scaling=2^") + to_string(scaling_bits)) {
            return ckks::encrypt(ckks::simd_encode(data, params), sk);
        };
    }

    for (auto [LOGN, scaling_bits] :
         {pair{12, 36}, pair{13, 43}, pair{14, 48}, pair{15, 55}}) {
        u64 N = 1 << LOGN;
        auto params = ckks::create_params(N, scaling_bits);
        cache_ntt_factors_strict(LOGN, params.moduli);

        CkksSk sk(params);
        std::vector<cc_double> data(N / 2);
        auto ct = ckks::encrypt(ckks::simd_encode(data, params), sk);
        BENCHMARK(string("CKKS decrypt+decode / N=") + to_string(N) +
                  string(" / scaling=2^") + to_string(scaling_bits)) {
            return ckks::simd_decode(ckks::decrypt(ct, sk));
        };
    }
}
