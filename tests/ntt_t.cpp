#include "catch2/catch.hpp"
#include "common/ntt.h"
#include "common/rnspolynomial.h"
#include <random>

using namespace hehub;
using SimplePoly = RnsPolynomial::ComponentData;

u64 powmod(u64 modulus, u64 base, size_t index) {
    u64 power = 1;
    size_t mask = 1;
    while (mask <= index) {
        mask <<= 1;
    }
    mask >>= 1;
    while (mask) {
        power = (u128)power * power % modulus;
        if (mask & index) {
            power = (u128)power * base % modulus;
        }
        mask >>= 1;
    }
    return power;
}

u64 bit_rev(u64 x, const size_t bit_len) {
    u64 mh = (1 << (bit_len - 1));
    u64 ml = 1;
    for (; mh >= ml; mh >>= 1, ml <<= 1) {
        if ((mh & x) && (ml & x)) {
        } else if (!(mh & x) && !(ml & x)) {
        } else if (!(mh & x) && (ml & x)) {
            x |= mh;
            x ^= ml;
        } else {
            x |= ml;
            x ^= mh;
        }
    }
    return x;
}

TEST_CASE("ntt", "[.]") {
    const size_t LOGN = GENERATE(11, 13, 15);
    const size_t N = 1 << LOGN;
    const auto [Q, R] = 
        GENERATE(std::make_pair<u64, u64>(260898817ULL, 5),
                 std::make_pair<u64, u64>(35184358850561ULL, 3), 
                 std::make_pair<u64, u64>(36028796997599233ULL, 5),
                 std::make_pair<u64, u64>(576460752272228353ULL, 5));

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

        for (size_t i = 0; i < N; i++) {
            REQUIRE(poly[i] < 2 * Q);
            REQUIRE(poly[i] % Q == 1);
        }
    }
    SECTION("just x") {
        for (int i = 0; i < N; i++) {
            poly[i] = 0;
        }
        poly[1] = 1;

        ntt_negacyclic_inplace_lazy(LOGN, Q, poly.data());        
        
        for (size_t i = 0; i < N; i++) {
            REQUIRE(poly[i] < 2 * Q);
        }

        auto root = powmod(Q, R, (Q - 1) / N / 2);
        for (size_t i = 0; i < N; i = i * 13 + 1) {
            auto curr_index = 2 * bit_rev(i, LOGN) + 1;
            auto curr_power = powmod(Q, root, curr_index);
            REQUIRE(poly[i] % Q == curr_power);
        }
    }
    SECTION("random") {
        for (int i = 0; i < N; i++) {
            poly[i] = distribution(generator);
        }
        auto poly_copy(poly);

        ntt_negacyclic_inplace_lazy(LOGN, Q, poly.data());

        for (size_t i = 0; i < N; i++) {
            REQUIRE(poly[i] < 2 * Q);
        }

        // naive NTT on selected points
        auto root = powmod(Q, R, (Q - 1) / N / 2);
        for (size_t i = 0; i < N; i = i * 13 + 1) {
            u64 curr_x = powmod(Q, root, 2 * bit_rev(i, LOGN) + 1);
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

TEST_CASE("ntt and intt") {
    auto LOGN = GENERATE(4, 7, 13, 14, 15);
    u64 N = 1 << LOGN;
    u64 Q = GENERATE(260898817ULL, 35184358850561ULL, 36028796997599233ULL,
                     576460752272228353ULL);

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

        for (int i = 0; i < N; i++) {
            REQUIRE(poly[i] < 2 * Q);
            if (poly[i] >= Q) {
                poly[i] -= Q;
            }
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

        for (int i = 0; i < N; i++) {
            REQUIRE(poly[i] < 2 * Q);
            if (poly[i] >= Q) {
                poly[i] -= Q;
            }
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

        for (int i = 0; i < N; i++) {
            REQUIRE(poly[i] < 2 * Q);
            if (poly[i] >= Q) {
                poly[i] -= Q;
            }
        }

        REQUIRE(poly == poly_copy);
    }
}
