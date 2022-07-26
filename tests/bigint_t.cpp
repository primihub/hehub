#include "catch2/catch.hpp"
#include "common/bigint.h"
#include "common/bigintpoly.h"
#include <iostream>
#include <sstream>

using namespace hehub;

TEST_CASE("big int") {
    UBInt first("123456789");
    UBInt second(123456789);
    REQUIRE(first == second);

    UBInt third("6666666666666666666666");
    UBInt fourth("8888888888888888888888");
    REQUIRE(third < fourth);
    UBInt fifth("12345678910111213141516171819");
    REQUIRE(fifth > fourth);

    // Printing all the numbers
    std::ostringstream oss;
    oss << first;
    REQUIRE(oss.str() == "123456789");
    oss.str("");
    oss << second;
    REQUIRE(oss.str() == "123456789");
    oss.str("");
    oss << third;
    REQUIRE(oss.str() == "6666666666666666666666");
    oss.str("");
    oss << fourth;
    REQUIRE(oss.str() == "8888888888888888888888");
    oss.str("");
    oss << fifth;
    REQUIRE(oss.str() == "12345678910111213141516171819");

    // Arithmetics
    first++;
    oss.str("");
    oss << first;
    REQUIRE(oss.str() == "123456790");
    UBInt sum = third + fourth;
    oss.str("");
    oss << sum;
    REQUIRE(oss.str() == "15555555555555555555554");
    UBInt product;
    product = second * third;
    oss.str("");
    oss << product;
    REQUIRE(oss.str() == "823045259999999999999917695474");
}

TEST_CASE("CRT composition") {
    std::vector<u64> moduli;
    std::vector<u64> remainders;
    SECTION("Single modulus") {
        moduli = std::vector<u64>{7};
        remainders = std::vector<u64>{3};
    }
    SECTION("multiple moduli") {
        moduli = std::vector<u64>{3, 5, 7};
        remainders = std::vector<u64>{2, 3, 4};
    }
    SECTION("big moduli") {
        moduli =
            std::vector<u64>{0x3ffffffffffe5, 0x3ffffffffffdd, 0x3ffffffffffcd};
        remainders = std::vector<u64>{123456789, 666666666, 888888888};
    }
    CRTComposer crt_composer(moduli);
    UBInt composed = crt_composer.compose(remainders);
    for (size_t i = 0; i < moduli.size(); i++) {
        REQUIRE(composed % UBInt(moduli[i]) == UBInt(remainders[i]));
    }
}

TEST_CASE("big int poly") {
    // Generate a pseudo-random RNS polynomial.
    auto moduli_vec =
        std::vector<u64>{0x3ffffffffffe5, 0x3ffffffffffdd, 0x3ffffffffffcd};
    RnsPolynomial rns_poly(3, 3, moduli_vec);
    auto seed = 42;
    auto poly_len = rns_poly.poly_len();
    auto components = rns_poly.component_count();
    for (size_t j = 0; j < components; j++) {
        auto curr_modulus = moduli_vec[j];
        for (size_t i = 0; i < poly_len; i++) {
            seed ^= seed * 498672493528838 + 93875838578748;
            rns_poly[j][i] = seed % curr_modulus;
        }
    }

    // Compose the RNS polynomial into a big-int polynomial, and check the
    // correctness.
    UBigIntPoly big_int_poly(rns_poly);
    REQUIRE(big_int_poly.poly_len() == poly_len);
    for (int j = 0; j < components; j++) {
        auto curr_modulus_big = UBInt(moduli_vec[j]);
        for (int i = 0; i < poly_len; i++) {
            auto curr_remainder_coeff_big = UBInt(rns_poly[j][i]);
            REQUIRE(big_int_poly[i] % curr_modulus_big ==
                    curr_remainder_coeff_big);
        }
    }
}
