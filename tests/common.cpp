#include "catch2/catch.hpp"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include "common/permutation.h"
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

TEST_CASE("bit rev", "[.]") {
    REQUIRE(__bit_rev_naive_16(12345, 14) == __bit_rev_naive(12345, 14));
    REQUIRE(__bit_rev_naive_16(12345, 15) == __bit_rev_naive(12345, 15));
    REQUIRE(__bit_rev_naive_16(12345, 16) == __bit_rev_naive(12345, 16));

#ifdef HEHUB_DEBUG
    REQUIRE_NOTHROW(__bit_rev_naive(12345, 64));
    REQUIRE_THROWS(__bit_rev_naive(12345, -1));
    REQUIRE_THROWS(__bit_rev_naive(12345, 10000000));
    REQUIRE_THROWS(__bit_rev_naive(12345, 13));

    REQUIRE_NOTHROW(__bit_rev_naive_16(12345, 16));
    REQUIRE_THROWS(__bit_rev_naive_16(12345, -1));
    REQUIRE_THROWS(__bit_rev_naive_16(12345, 10000000));
    REQUIRE_THROWS(__bit_rev_naive_16(12345, 13));
#endif
}

/// Infinity norm on a simple (not RNS) polynomial
u64 simple_inf_norm(const RnsPolynomial &poly) {
    if (poly.component_count() != 1) {
        throw std::invalid_argument("poly");
    }
    if (poly.rep_form != PolyRepForm::coeff) {
        throw std::invalid_argument("poly");
    }

    u64 q = poly.modulus_at(0);
    u64 half_q = q / 2;
    u64 norm = 0;
    for (const auto &coeff : poly[0]) {
        if (coeff >= q) {
            throw std::logic_error("Not reduced strictly.");
        }
        if (coeff < half_q) {
            norm = std::max(coeff, norm);
        } else {
            norm = std::max(q - coeff, norm);
        }
    }

    return norm;
}

TEST_CASE("automorphism") {
    SECTION("tables") {
        size_t logn = 3;
        size_t n = 1 << logn;
        auto dlog_table = dlog_mod_2power_table(logn);
        std::vector<size_t> dlog_table_compact(n);
        for (size_t i = 0; i < n; i++) {
            dlog_table_compact[i] = dlog_table[2 * i + 1];
        }

        REQUIRE(dlog_table_compact ==
                std::vector<size_t>{0, 1, 4, 5, 2, 3, 6, 7});
    }
    SECTION("involution") {
        u64 q = 65537;
        RnsPolynomial poly(4096, 1, {q});

        // random polynomial with not too large norm
        u64 seed = 42;
        for (auto &v : poly[0]) {
            seed = seed * 4985348 + 93479384;
            v = seed % (q / 10);
        }
        ntt_negacyclic_inplace_lazy(poly);

        // check the involution
        auto involuted = involute(poly);
        REQUIRE(involute(involuted) == poly);

        // check the boundness property
        intt_negacyclic_inplace(poly);
        intt_negacyclic_inplace(involuted);
        REQUIRE(simple_inf_norm(poly) == simple_inf_norm(involuted));
    }
    SECTION("cycles") {
        u64 q = 65537;
        size_t poly_len = 4096;
        RnsPolynomial poly(4096, 1, {q});

        // random polynomial with not too large norm
        u64 seed = 42;
        for (auto &v : poly[0]) {
            seed = seed * 4985348 + 93479384;
            v = seed % (q / 10);
        }
        ntt_negacyclic_inplace_lazy(poly);

        // check the cycles
        auto one_step = cycle(poly, 1);
        auto two_step = cycle(poly, 2);
        REQUIRE(cycle(one_step, poly_len / 2 - 1) == poly);
        REQUIRE(cycle(one_step, 1) == two_step);

        // check the boundness property
        intt_negacyclic_inplace(poly);
        intt_negacyclic_inplace(one_step);
        intt_negacyclic_inplace(two_step);
        REQUIRE(simple_inf_norm(poly) == simple_inf_norm(one_step));
        REQUIRE(simple_inf_norm(poly) == simple_inf_norm(two_step));
    }
}
