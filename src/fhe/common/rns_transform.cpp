#include "rns_transform.h"
#include "bigintpoly.h"
#include "mod_arith.h"
#include "range/v3/view/zip.hpp"
#include <numeric>

using namespace ranges::views;

namespace hehub {

RnsPolynomial
rns_base_transform_from_single(const RnsPolynomial &input_rns_poly,
                               const std::vector<u64> &new_moduli) {
    auto old_modulus = input_rns_poly.modulus_at(0);
    auto half_old_modulus = old_modulus / 2;
    auto dimension = input_rns_poly.dimension();

    RnsPolyParams output_dim{dimension, new_moduli.size(), new_moduli};
    RnsPolynomial result(output_dim);
    auto &input_poly = input_rns_poly[0];
    for (auto [component, modulus] : zip(result, new_moduli)) {
        auto modulus_multiple = (old_modulus / modulus + 1) * modulus;
        for (auto [component_coeff, input_coeff] : zip(component, input_poly)) {
            if (input_coeff < half_old_modulus) {
                component_coeff = input_coeff;
            } else {
                component_coeff = modulus_multiple - old_modulus + input_coeff;
            }
        }

        if (modulus < old_modulus) {
            batched_barrett_lazy(modulus, dimension, component.data());
        }
    }

    return result;
}

RnsPolynomial rns_base_transform_to_single(const RnsPolynomial &input_rns_poly,
                                           const u64 new_modulus) {
    auto dimension = input_rns_poly.dimension();
    auto &old_moduli = input_rns_poly.modulus_vec();
    RnsPolyParams params{dimension, 1, std::vector{new_modulus}};
    RnsPolynomial result(params);
    auto &result_poly = result[0];

    // Check if the coefficients are smaller than each old modulus. If so then
    // their composed form is equivalent to any one component.
    bool small_coeff = true;
    auto first_old_mod = old_moduli[0];
    u64 half_first_old_mod = first_old_mod / 2;
    for (size_t i = 0; i < dimension; i++) {
        for (size_t k = 1; k < input_rns_poly.component_count(); k++) {
            if (input_rns_poly[0][i] < half_first_old_mod &&
                input_rns_poly[k][i] != input_rns_poly[0][i]) {
                small_coeff = false;
                break;
            }
            if (input_rns_poly[0][i] >= half_first_old_mod &&
                old_moduli[k] - input_rns_poly[k][i] !=
                    old_moduli[0] - input_rns_poly[0][i]) {
                small_coeff = false;
                break;
            }
        }
        if (!small_coeff) {
            break;
        }
    }
    if (small_coeff) {
        u64 new_modulus_multiple =
            (first_old_mod / new_modulus + 1) * new_modulus;
        for (auto [result_coeff, input_coeff] :
             zip(result_poly, input_rns_poly[0])) {
            if (input_coeff < half_first_old_mod) {
                result_coeff = input_coeff;
            } else {
                result_coeff =
                    new_modulus_multiple - first_old_mod + input_coeff;
            }
        }
        batched_barrett(new_modulus, dimension, result_poly.data());
        return result;
    }

    // If the coefficients are not all smaller than each modulus, then we need
    // to perform a CRT composition.
    UBigIntPoly big_int_poly(input_rns_poly);
    auto big_int_new_modulus = UBInt(new_modulus);
    auto big_int_old_modulus = UBInt(1);
    for (const auto &mod : input_rns_poly.modulus_vec()) {
        big_int_old_modulus *= UBInt(mod);
    }
    auto half_big_int_old_modulus = big_int_old_modulus / UBInt(2);
    for (size_t i = 0; i < dimension; i++) {
        if (big_int_poly[i] < half_big_int_old_modulus) {
            result_poly[i] = to_u64(big_int_poly[i] % big_int_new_modulus);
        } else {
            auto abs_val = big_int_old_modulus - big_int_poly[i];
            abs_val %= big_int_new_modulus;
            result_poly[i] = new_modulus - to_u64(abs_val);
        }
    }
    return result;
}

RnsPolynomial rns_base_transform(RnsPolynomial input_rns_poly,
                                 const std::vector<u64> &new_moduli) {
    if (input_rns_poly.rep_form == PolyRepForm::value) {
        throw std::logic_error("Trying to perform RNS base transformation "
                               "on NTT values.");
    }

    // Reduce the coefficients strictly to avoid errors caused by redundant
    // multiples of original modulus
    reduce_strict(input_rns_poly);

    if (input_rns_poly.component_count() == 1) {
        return rns_base_transform_from_single(input_rns_poly, new_moduli);
    } else if (new_moduli.size() == 1) {
        return rns_base_transform_to_single(input_rns_poly, new_moduli[0]);
    } else {
        throw "under development";
    }
}

} // namespace hehub
