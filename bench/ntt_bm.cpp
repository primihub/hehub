// #include "catch2/catch.hpp"
// #include "fhe/common/ntt.h"
// #include "fhe/common/rns.h"
// #include <random>

// using namespace hehub;

// TEST_CASE("benchmark ntt and intt") {
//     u64 Q = 576460752272228353ULL;

//     std::random_device rd;
//     std::default_random_engine generator(rd());
//     std::uniform_int_distribution<u64> distribution(0, Q - 1);

//     for (auto LOGN: {10, 11, 12, 13, 14, 15}) {
//         u64 N = 1 << LOGN;
//         cache_ntt_factors_strict(LOGN, std::vector{Q});

//         u64 poly[N];
//         for (int i = 0; i < N; i++) {
//             poly[i] = distribution(generator);
//         }
        
//         BENCHMARK(std::string("NTT / len=") + std::to_string(N)) {
//             return ntt_negacyclic_inplace_lazy(LOGN, Q, poly);
//         };
//     }

//     for (auto LOGN: {10, 11, 12, 13, 14, 15}) {
//         u64 N = 1 << LOGN;

//         u64 poly[N];
//         for (int i = 0; i < N; i++) {
//             poly[i] = distribution(generator);
//         }
        
//         BENCHMARK(std::string("INTT / len=") + std::to_string(N)) {
//             return intt_negacyclic_inplace_lazy(LOGN, Q, poly);
//         };
//     }
// }
