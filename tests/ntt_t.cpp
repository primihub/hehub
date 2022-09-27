#include "catch2/catch.hpp"
#include "common/mod_arith.h"
#include "common/ntt.h"
#include "common/permutation.h"
#include "common/rnspolynomial.h"
#include <numeric>
#include <random>

namespace hehub {
u64 __pow_mod(u64 modulus, u64 base, size_t index);

u64 __get_2nth_unity_root(u64 modulus, u64 n);
} // namespace hehub

using namespace hehub;
using SimplePoly = RnsPolynomial::ComponentData;

TEST_CASE("ntt", "[.]") {
    const size_t LOGN = GENERATE(11, 13, 15);
    const size_t N = 1 << LOGN;
    const u64 Q = GENERATE(260898817ULL, 35184358850561ULL,
                           36028796997599233ULL, 576460752272228353ULL);

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_int_distribution<u64> distribution(0, Q - 1);

    SimplePoly poly(N);

    SECTION("one") {
        for (int i = 0; i < N; i++) {
            poly[i] = 0;
        }
        poly[0] = 1;

        ntt_negacyclic_inplace_lazy(LOGN, Q, poly.data());

        for (const auto &value : poly) {
            REQUIRE(value < 2 * Q);
            REQUIRE(value % Q == 1);
        }
    }
    SECTION("just x") {
        for (auto &coeff : poly) {
            coeff = 0;
        }
        poly[1] = 1;

        ntt_negacyclic_inplace_lazy(LOGN, Q, poly.data());

        for (const auto &value : poly) {
            REQUIRE(value < 2 * Q);
        }

        auto root = __get_2nth_unity_root(Q, N);
        for (size_t i = 0; i < N; i = i * 13 + 1) {
            auto curr_index = 2 * __bit_rev_naive(i, LOGN) + 1;
            auto curr_power = __pow_mod(Q, root, curr_index);
            REQUIRE(poly[i] % Q == curr_power);
        }
    }
    SECTION("random") {
        for (auto &coeff : poly) {
            coeff = distribution(generator);
        }
        auto poly_copy(poly);

        ntt_negacyclic_inplace_lazy(LOGN, Q, poly.data());

        for (const auto &value : poly) {
            REQUIRE(value < 2 * Q);
        }

        // naive NTT on selected points
        auto root = __get_2nth_unity_root(Q, N);
        for (size_t i = 0; i < N; i = i * 13 + 1) {
            u64 curr_x = __pow_mod(Q, root, 2 * __bit_rev_naive(i, LOGN) + 1);
            u64 acc = 0;
            u64 x_power = 1;
            for (int j = 0; j < N; j++) {
                acc += (u128)poly_copy[j] * x_power % Q;
                acc %= Q;
                x_power = (u128)x_power * curr_x % Q;
            }

            REQUIRE(poly[i] % Q == acc);
        }
    }
}

TEST_CASE("ntt round trip") {
    auto LOGN = GENERATE(4, 7, 13, 14, 15);
    u64 N = 1 << LOGN;
    u64 Q = GENERATE(65537ULL, 260898817ULL, 35184358850561ULL,
                     36028796997599233ULL, 576460752272228353ULL);

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_int_distribution<u64> distribution(0, Q - 1);

    SimplePoly poly(N);

    SECTION("one") {
        for (int i = 0; i < N; i++) {
            poly[i] = 0;
        }
        poly[0] = 1;

        SimplePoly poly_copy = poly;

        ntt_negacyclic_inplace_lazy(LOGN, Q, poly.data());
        intt_negacyclic_inplace_lazy(LOGN, Q, poly.data());

        auto check_range = [=](bool status, u64 coeff) {
            return status && (coeff < 2 * Q);
        };
        REQUIRE(
            std::accumulate(poly.data(), poly.data() + N, true, check_range));
        for (int i = 0; i < N; i++) {
            poly[i] -= (poly[i] >= Q) ? Q : 0;
        }

        REQUIRE(poly == poly_copy);
    }
    SECTION("just x") {
        for (int i = 0; i < N; i++) {
            poly[i] = 0;
        }
        poly[1] = 1;

        SimplePoly poly_copy = poly;

        ntt_negacyclic_inplace_lazy(LOGN, Q, poly.data());
        intt_negacyclic_inplace_lazy(LOGN, Q, poly.data());

        auto check_range = [=](bool status, u64 coeff) {
            return status && (coeff < 2 * Q);
        };
        REQUIRE(
            std::accumulate(poly.data(), poly.data() + N, true, check_range));
        for (int i = 0; i < N; i++) {
            poly[i] -= (poly[i] >= Q) ? Q : 0;
        }

        REQUIRE(poly == poly_copy);
    }
    SECTION("random") {
        for (int i = 0; i < N; i++) {
            poly[i] = distribution(generator);
        }

        SimplePoly poly_copy = poly;

        ntt_negacyclic_inplace_lazy(LOGN, Q, poly.data());
        intt_negacyclic_inplace_lazy(LOGN, Q, poly.data());

        auto check_range = [=](bool status, u64 coeff) {
            return status && (coeff < 2 * Q);
        };
        REQUIRE(
            std::accumulate(poly.data(), poly.data() + N, true, check_range));
        for (int i = 0; i < N; i++) {
            poly[i] -= (poly[i] >= Q) ? Q : 0;
        }

        REQUIRE(poly == poly_copy);
    }
    SECTION("on rns poly") {
        RnsPolynomial rns_poly(N, 1, std::vector{Q});
        for (int i = 0; i < N; i++) {
            rns_poly[0][i] = distribution(generator);
        }
        auto rns_poly_copy(rns_poly);

        ntt_negacyclic_inplace_lazy(rns_poly);
        intt_negacyclic_inplace_lazy(rns_poly);
        strict_reduce(rns_poly);

        REQUIRE(rns_poly == rns_poly_copy);
    }
}
