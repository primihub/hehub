#include "sampling.h"
#include "ntt.h"
#include "range/v3/view/zip.hpp"
#include "rns.h"
#include "type_defs.h"
#include <random>

using namespace ranges::views;

namespace hehub {

std::random_device rand_dvc;
std::default_random_engine rand_engine;
std::uniform_int_distribution rand_ternary((i8)-1, (i8)1);

RnsPolynomial get_rand_ternary_poly(const RlweParams &params) {
    RnsPolynomial tern_poly(params);
    auto dimension = params.dimension;

    // Sampling.
    std::vector<i8> ternary_integers(dimension);
    for (auto &t : ternary_integers) {
        t = rand_ternary(rand_engine);
    }

    // Transform to RNS representation.
    for (auto [component, modulus] : zip(tern_poly, tern_poly.modulus_vec())) {
        for (auto [coeff, ternary_int] : zip(component, ternary_integers)) {
            coeff = modulus + (u64)ternary_int;
            coeff -= (coeff >= modulus) ? modulus : 0;
        }
    }

    // Transform to NTT values.
    ntt_negacyclic_inplace_lazy(tern_poly);
    return tern_poly;
}

RnsPolynomial get_rand_uniform_poly(const RlweParams &params,
                                    PolyRepForm form) {
    auto dimension = params.dimension;
    RnsPolynomial rand_rns_poly(params);

    // Sampling.
    for (auto [component, modulus] :
         zip(rand_rns_poly, rand_rns_poly.modulus_vec())) {
        std::uniform_int_distribution uni_mod((u64)0, modulus - 1);
        for (auto &coeff : component) {
            // here "coeff" can also mean NTT value
            coeff = uni_mod(rand_engine);
        }
    }

    // Set the representation form label.
    rand_rns_poly.rep_form = form;

    return rand_rns_poly;
}

RnsPolynomial get_rand_gaussian_poly(const RlweParams &params,
                                     double std_dev) {
    auto dimension = params.dimension;
    RnsPolynomial gaussian_poly(params);

    auto bound = std_dev * 6;
    std::normal_distribution<double> rand_gaussian(0, std_dev);

    // Sampling.
    std::vector<double> gassians(dimension);
    for (auto &g : gassians) {
        do {
            g = rand_gaussian(rand_engine);
        } while (std::abs(g) > bound);
    }

    // Transform to RNS representation.
    for (auto [component, modulus] :
         zip(gaussian_poly, gaussian_poly.modulus_vec())) {
        for (auto [coeff, gaussian] : zip(component, gassians)) {
            coeff = modulus + gaussian;
            coeff -= (coeff >= modulus) ? modulus : 0;
        }
    }

    // Transform to NTT values.
    ntt_negacyclic_inplace_lazy(gaussian_poly);
    return gaussian_poly;
}

RnsPolynomial get_zero_poly(const RlweParams &params, PolyRepForm form) {
    RnsPolynomial rns_poly(params);
    rns_poly.rep_form = form;
    for (auto &component : rns_poly) {
        std::fill(component.begin(), component.end(), 0);
    }
    return rns_poly;
}

} // namespace hehub
