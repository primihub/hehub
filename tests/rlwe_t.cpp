#include "catch2/catch.hpp"
#include "primitives/rlwe.h"
#include <numeric>

using namespace hehub;

TEST_CASE("rlwe") {
    std::vector<u64> moduli{35184358850561ULL, 36028796997599233ULL,
                            576460752272228353ULL};
    size_t poly_len = 4096;
    PolyDimensions poly_dim{poly_len, moduli.size(), moduli};
    RlwePt pt(poly_dim);
    RlweSk sk(poly_dim);

    // plaintext data
    const i64 DATUM_TEST = 1000000000;
    for (auto &component_poly : pt) {
        for (auto &datum : component_poly) {
            datum = DATUM_TEST;
        }
    }

    // encrypt & decrypt
    RlweCt ct = encrypt(pt, sk);
    RlwePt pt_recovered = decrypt(ct, sk);

    // check
    auto component_count = pt.component_count();
    REQUIRE(pt_recovered.component_count() == component_count);
    auto check_if_close = [=](auto status, auto coeff) {
        return status && (std::abs(i64(coeff) - DATUM_TEST) < 20);
    };
    for (auto &component_poly : pt_recovered) {
        REQUIRE(std::accumulate(component_poly.begin(), component_poly.end(),
                                true, check_if_close));
    }
}
