// #define CATCH_CONFIG_MAIN
// #include "catch2/catch.hpp"

#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"
#include "fhe/common/ntt.h"
#include "fhe/ckks/ckks.h"

using namespace hehub;
using namespace std;

int main() {
    // double d = 1.0;
    // ankerl::nanobench::Bench().run("some double ops", [&] {
    //     d += 1.0 / d;
    //     if (d > 5.0) {
    //         d -= 5.0;
    //     }
    //     ankerl::nanobench::doNotOptimizeAway(d);
    // });
    for (auto [LOGN, scaling_bits] :
         {pair{12, 36}, pair{13, 43}, pair{14, 48}, pair{15, 55}}) {
        u64 N = 1 << LOGN;
        auto params = ckks::create_params(N, scaling_bits);
        cache_ntt_factors_strict(LOGN, params.moduli);

        CkksSk sk(params);
        auto rot_key = get_rot_key(sk, params.additional_mod, 1);
        std::vector<cc_double> data(N / 2);
        auto ct = ckks::encrypt(ckks::simd_encode(data, params), sk);
        ankerl::nanobench::Bench().timeUnit(1ms, "ms")
            .run(string("CKKS rotation / N=") + to_string(N) +
                 string(" / scaling=2^") + to_string(scaling_bits), [&] {
            auto ct_rotated = ckks::rotate(ct, rot_key, 1);
            ankerl::nanobench::doNotOptimizeAway(ct_rotated);
        });
    }
}
