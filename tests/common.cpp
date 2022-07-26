#include "catch2/catch.hpp"
#include "common/rnspolynomial.h"

using namespace hehub;

TEST_CASE("RNS poly") {
    RnsPolynomial r1(3, 12, std::vector<u64>{3, 5, 7});
    RnsPolynomial r2(3, 12, std::vector<u64>{3, 5, 7});
    RnsPolynomial r3(3, 12, std::vector<u64>{3, 5, 7});

    RnsPolynomial r4(r2);
    RnsPolynomial r5(std::move(r1));
    r3 = r2;
    r4 = std::move(r2);

    r3.add_components(std::vector<u64>{11});
    r4.remove_components();

    REQUIRE(r3.component_count() == 4);
    REQUIRE(r4.component_count() == 2);
    REQUIRE(r5.component_count() == 3);
}
