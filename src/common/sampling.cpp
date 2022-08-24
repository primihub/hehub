#include "sampling.h"
#include "ntt.h"
#include "rnspolynomial.h"
#include "type_defs.h"
#include <iostream>
#include <random>

namespace hehub {

std::random_device rand_dvc;
std::default_random_engine rand_engine;
std::uniform_int_distribution rand_ternary((i8)-1, (i8)1);

RnsPolynomial get_rand_ternary_poly(const PolyDimensions &poly_dim) {
    RnsPolynomial tern_poly(poly_dim);
    auto poly_len = poly_dim.poly_len;

    // Sampling.
    i8 temp[poly_len];
    for (auto &temp_item : temp) {
        temp_item = rand_ternary(rand_engine);
    }

    // Transform to RNS representation.
    auto mod_ptr = tern_poly.moduli_vec().begin();
    for (auto &component_poly : tern_poly) {
        auto curr_mod = *(mod_ptr++);
        auto temp_item_ptr = temp;
        for (auto &coeff : component_poly) {
            auto temp_item = u64(*(temp_item_ptr++));
            coeff = curr_mod + temp_item;
            coeff -= (coeff >= curr_mod) ? curr_mod : 0;
        }
    }

    // Transform to NTT values.
    ntt_negacyclic_inplace_lazy(tern_poly);
    return tern_poly;
}

RnsPolynomial get_rand_uniform_poly(const PolyDimensions &poly_dim,
                                    PolyRepForm form) {
    auto poly_len = poly_dim.poly_len;
    auto mod_ptr = poly_dim.moduli.begin();
    RnsPolynomial rand_rns_poly(poly_dim);

    // Sampling.
    for (auto &component_poly : rand_rns_poly) {
        std::uniform_int_distribution uni_mod((u64)0, *(mod_ptr++) - 1);
        for (auto &coeff : component_poly) {
            // here "coeff" can also be NTT value
            coeff = uni_mod(rand_engine);
        }
    }

    // Set the representation form label.
    rand_rns_poly.rep_form = form;

    return rand_rns_poly;
}

RnsPolynomial get_rand_gaussian_poly(const PolyDimensions &poly_dim,
                                     double std_dev) {
    auto mod_ptr = poly_dim.moduli.begin();
    auto poly_len = poly_dim.poly_len;
    RnsPolynomial gaussian_poly(poly_dim);

    auto bound = std_dev * 6;
    std::normal_distribution<double> rand_gaussian(0, std_dev);

    // Sampling.
    double temp[poly_len];
    for (auto &temp_item : temp) {
        do {
            temp_item = rand_gaussian(rand_engine);
        } while (std::abs(temp_item) > bound);
    }

    // Transform to RNS representation.
    for (auto &component_poly : gaussian_poly) {
        auto curr_mod = *(mod_ptr++);
        auto temp_item_ptr = temp;
        for (auto &coeff : component_poly) {
            auto temp_item = u64(*(temp_item_ptr++));
            coeff = curr_mod + temp_item;
            coeff -= (coeff >= curr_mod) ? curr_mod : 0;
        }
    }

    // Transform to NTT values.
    ntt_negacyclic_inplace_lazy(gaussian_poly);
    return gaussian_poly;
}

} // namespace hehub
