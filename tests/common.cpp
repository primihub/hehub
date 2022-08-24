#include "catch2/catch.hpp"
#include "common/rnspolynomial.h"

using namespace hehub;

TEST_CASE("RNS polynomial") {
    RnsPolynomial r1(4096, 3, std::vector<u64>{3, 5, 7});

    PolyDimensions poly_dim{4096, 3, std::vector<u64>{3, 5, 7}};
    RnsPolynomial r2(poly_dim);
    RnsPolynomial r3(poly_dim);

    RnsPolynomial r4(r2);
    RnsPolynomial r5(std::move(r1));
    r3 = r2;
    r4 = std::move(r2);

    r3.add_components(std::vector<u64>{11});
    r4.remove_components();

    REQUIRE(r3.component_count() == 4);
    REQUIRE(r4.component_count() == 2);
    REQUIRE(r5.component_count() == 3);

    REQUIRE_THROWS(RnsPolynomial(PolyDimensions{4096, 4, std::vector<u64>(3)}));
    REQUIRE_THROWS(RnsPolynomial(PolyDimensions{4095, 3, std::vector<u64>(3)}));
    REQUIRE_THROWS(RnsPolynomial(PolyDimensions{4097, 3, std::vector<u64>(3)}));
}
