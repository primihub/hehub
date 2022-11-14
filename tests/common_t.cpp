#include "catch2/catch.hpp"
#include "cereal/archives/binary.hpp"
#include "fhe/ckks/ckks.h"
#include "fhe/common/mod_arith.h"
#include "fhe/common/ntt.h"
#include "fhe/common/permutation.h"
#include "fhe/common/rns.h"
#include "fhe/common/sampling.h"
#include "fhe/primitives/rlwe.h"
#include <fstream>

using namespace hehub;

TEST_CASE("RNS polynomial") {
    RnsPolynomial r1(4096, 3, std::vector<u64>{3, 5, 7});

    RnsPolyParams params{4096, 3, std::vector<u64>{3, 5, 7}};
    RnsPolynomial r2(params);
    RnsPolynomial r3(params);

    RnsPolynomial r4(r2);
    RnsPolynomial r5(std::move(r1));
    r3 = r2;
    r4 = std::move(r2);

    r3.add_components(std::vector<u64>{11});
    r4.remove_components();

    REQUIRE(r3.component_count() == 4);
    REQUIRE(r4.component_count() == 2);
    REQUIRE(r5.component_count() == 3);

    REQUIRE_THROWS(RnsPolynomial(RnsPolyParams{4096, 4, std::vector<u64>(3)}));
    REQUIRE_THROWS(RnsPolynomial(RnsPolyParams{4095, 3, std::vector<u64>(3)}));
    REQUIRE_THROWS(RnsPolynomial(RnsPolyParams{4097, 3, std::vector<u64>(3)}));
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

TEST_CASE("sampling", "[.]") {
    u64 mod = 65537;
    RnsPolyParams params{4096, 1, std::vector<u64>{mod}};

    SECTION("ternary") {
        auto tern_poly = get_rand_ternary_poly(params);
        intt_negacyclic_inplace(tern_poly);

        for (const auto &coeff : tern_poly[0]) {
            REQUIRE(((coeff == 0) || (coeff == 1) || (coeff == mod - 1)));
        }

        bool all_zero = true;
        for (const auto &coeff : tern_poly[0]) {
            if (coeff != 0) {
                all_zero = false;
            }
        }
        REQUIRE_FALSE(all_zero);
    }
    SECTION("gaussian") {
        auto gauss_poly = get_rand_gaussian_poly(params);
        intt_negacyclic_inplace(gauss_poly);

        for (const auto &coeff : gauss_poly[0]) {
            REQUIRE(((coeff < 20) || (coeff > mod - 20)));
        }

        bool all_zero = true;
        for (const auto &coeff : gauss_poly[0]) {
            if (coeff != 0) {
                all_zero = false;
            }
        }
        REQUIRE_FALSE(all_zero);
    }
    SECTION("uniform") {
        auto uniform_poly = get_rand_uniform_poly(params);
        intt_negacyclic_inplace(uniform_poly);

        bool all_zero = true;
        for (const auto &coeff : uniform_poly[0]) {
            if (coeff != 0) {
                all_zero = false;
            }
        }
        REQUIRE_FALSE(all_zero);
    }
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

template <typename T>
bool all_close(const std::vector<T> &vec1, const std::vector<T> &vec2,
               double eps) {
    for (size_t i = 0; i < vec1.size(); i++) {
        if (i >= vec2.size() || std::abs(vec1[i] - vec2[i]) > eps) {
            return false;
        }
    }
    return true;
};

TEST_CASE("automorphism") {
    SECTION("involution") {
        u64 q = 65537;
        size_t dimension = 8;
        RnsPolynomial poly(dimension, 1, {q});

        // random polynomial with not too large norm
        u64 seed = 42;
        for (auto &v : poly[0]) {
            seed = seed * 4985348 + 93479384;
            v = seed % (q / 10);
        }
        ntt_negacyclic_inplace_lazy(poly);

        // check the involution
        auto poly_involved = involution(poly);
        REQUIRE(involution(poly_involved) == poly);

        // check the boundness property
        intt_negacyclic_inplace(poly);
        intt_negacyclic_inplace(poly_involved);
        REQUIRE(simple_inf_norm(poly) == simple_inf_norm(poly_involved));
    }
    SECTION("cycles") {
        u64 q = 65537;
        size_t dimension = 8;
        RnsPolynomial poly(dimension, 1, {q});

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
        REQUIRE(cycle(one_step, dimension / 2 - 1) == poly);
        REQUIRE(cycle(one_step, 1) == two_step);

        // check the boundness property
        intt_negacyclic_inplace(poly);
        intt_negacyclic_inplace(one_step);
        intt_negacyclic_inplace(two_step);
        REQUIRE(simple_inf_norm(poly) == simple_inf_norm(one_step));
        REQUIRE(simple_inf_norm(poly) == simple_inf_norm(two_step));
    }
    SECTION("involution on plain") {
        size_t dimension = 8;
        CkksParams params = create_params(dimension, {55});
        params.initial_scaling_factor = pow(2.0, 50);
        auto data_count = dimension / 2;
        std::vector<cc_double> plain_data(data_count);
        std::vector<cc_double> data_conj;
        std::default_random_engine generator;
        std::normal_distribution<double> data_dist(0, 1);
        for (auto &d : plain_data) {
            d = {data_dist(generator), data_dist(generator)};
            data_conj.push_back(std::conj(d));
        }

        auto pt = ckks::simd_encode(plain_data, params);
        ntt_negacyclic_inplace_lazy(pt);
        CkksPt pt_involved = involution(pt);
        pt_involved.scaling_factor = pt.scaling_factor;
        intt_negacyclic_inplace(pt_involved);
        auto data_recovered = ckks::simd_decode<cc_double>(pt_involved);

        REQUIRE(all_close(data_recovered, data_conj, pow(2.0, -45)));
    }
    SECTION("cycle on plain") {
        size_t dimension = 8;
        CkksParams params = create_params(dimension, {55});
        params.initial_scaling_factor = pow(2.0, 50);
        auto data_count = dimension / 2;
        std::vector<cc_double> plain_data(data_count);
        std::vector<cc_double> data_rot(data_count);
        std::default_random_engine generator;
        std::normal_distribution<double> data_dist(0, 1);
        for (auto &d : plain_data) {
            d = {data_dist(generator), data_dist(generator)};
        }
        size_t step = GENERATE(1, 2, 3);
        for (size_t i = 0; i < data_count; i++) {
            data_rot[(i + step) % data_count] = plain_data[i];
        }

        auto pt = ckks::simd_encode(plain_data, params);
        ntt_negacyclic_inplace_lazy(pt);
        CkksPt pt_cycled = cycle(pt, step);
        pt_cycled.scaling_factor = pt.scaling_factor;
        intt_negacyclic_inplace(pt_cycled);
        auto data_recovered = ckks::simd_decode<cc_double>(pt_cycled);

        REQUIRE(all_close(data_recovered, data_rot, pow(2.0, -45)));
    }
}

TEST_CASE("serialization") {
    std::stringstream ss;
    cereal::BinaryOutputArchive oarchive(ss);
    cereal::BinaryInputArchive iarchive(ss);

    SECTION("rns poly component") {
        using SimplePoly = RnsPolynomial::ComponentData;

        SimplePoly simple_poly(8);
        simple_poly.save(oarchive);

        SimplePoly simple_poly2;
        simple_poly2.load(iarchive);

        REQUIRE(simple_poly2 == simple_poly);
    }
    SECTION("rns poly") {
        RlweParams params{4096, 3, {3, 5, 7}};
        RnsPolynomial rns_poly(params);
        rns_poly.rep_form = GENERATE(PolyRepForm::coeff, PolyRepForm::value);
        save(oarchive, rns_poly);

        RnsPolynomial rns_poly2;
        load(iarchive, rns_poly2);

        REQUIRE(rns_poly2 == rns_poly);
    }
    SECTION("rlwe params") {
        RlweParams params = create_params(8, {30, 30});
        save(oarchive, params);

        RlweParams params2;
        load(iarchive, params2);

        REQUIRE(params2 == params);
    }
#ifdef HEHUB_DEBUG_FHE
    SECTION("rlwe sk") {
        RlweParams params = create_params(8, {30, 30});
        RlweSk sk(params);
        save(oarchive, sk);

        RlweSk sk2;
        load(iarchive, sk2);

        REQUIRE(sk2 == sk);
    }
#endif
    SECTION("rlwe") {
        RlweParams params = create_params(8, {30, 30});
        RlweSk sk(params);
        RlweCt ct = get_rlwe_sample(sk);
        save(oarchive, ct);

        RlweCt ct2;
        load(iarchive, ct2);

        REQUIRE(ct2 == ct);
    }
    SECTION("ckks params") {
        CkksParams params = ckks::create_params(8, {30, 30}, 30, pow(2.0, 30));
        save(oarchive, params);

        CkksParams params2;
        load(iarchive, params2);

        REQUIRE(params2 == params);
    }
    SECTION("ckks ct") {
        CkksParams params = ckks::create_params(8, {30, 30}, 30, pow(2.0, 30));
        CkksSk sk(params);
        CkksCt ct = get_rlwe_sample(sk);
        ct.scaling_factor = params.initial_scaling_factor;
        save(oarchive, ct);

        CkksCt ct2;
        load(iarchive, ct2);

        REQUIRE(ct2 == ct);
    }
}
