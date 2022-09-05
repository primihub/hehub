#include "rns_transform.h"
#include "bigintpoly.h"
#include "mod_arith.h"
#include <cassert>
#include <numeric>

namespace hehub {

RnsPolynomial
rns_base_transform_from_single(const RnsPolynomial input_rns_poly,
                               const std::vector<u64> &new_moduli) {
    auto old_modulus = input_rns_poly.modulus_at(0);
    auto poly_len = input_rns_poly.poly_len();
    PolyDimensions poly_dim{poly_len, new_moduli.size(), new_moduli};
    RnsPolynomial result(poly_dim);

    auto mod_ptr = new_moduli.begin();
    for (auto &component_poly : result) {
        component_poly = input_rns_poly[0];

        auto component_mod = *(mod_ptr++);
        if (component_mod < old_modulus) {
            batched_barrett_lazy(component_mod, poly_len,
                                 component_poly.data());
        }
    }

    return result;
}

RnsPolynomial rns_base_transform_to_single(const RnsPolynomial input_rns_poly,
                                           const u64 new_modulus) {
    auto poly_len = input_rns_poly.poly_len();
    PolyDimensions poly_dim{poly_len, 1, std::vector{new_modulus}};
    RnsPolynomial result(poly_dim);
    auto &result_poly = result[0];

    for (size_t k = 0; k < input_rns_poly.component_count(); k++) {
        auto component_mod = input_rns_poly.modulus_at(k);
        if (component_mod > new_modulus) {
            u64 half_component_mod = component_mod / 2;
            u64 new_modulus_multiple =
                (component_mod / new_modulus + 1) * new_modulus;
            for (size_t i = 0; i < poly_len; i++) {
                if (input_rns_poly[k][i] < half_component_mod) {
                    result_poly[i] = input_rns_poly[k][i];
                } else {
                    result_poly[i] = new_modulus_multiple - component_mod +
                                     input_rns_poly[k][i];
                }
            }
            batched_barrett(new_modulus, poly_len, result_poly.data());
            return result;
        }
    }

    // If new moduli are all smaller than the old, then we need to perform a CRT
    // composition.
    UBigIntPoly big_int_poly(input_rns_poly);
    auto big_int_new_modulus = UBInt(new_modulus);
    auto big_int_old_modulus = UBInt(1);
    for (const auto &mod : input_rns_poly.moduli_vec()) {
        big_int_old_modulus *= UBInt(mod);
    }
    auto half_big_int_old_modulus = big_int_old_modulus / UBInt(2);
    for (size_t i = 0; i < poly_len; i++) {
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
    strict_reduce(input_rns_poly);

    if (input_rns_poly.component_count() == 1) {
        return rns_base_transform_from_single(input_rns_poly, new_moduli);
    } else if (new_moduli.size() == 1) {
        return rns_base_transform_to_single(input_rns_poly, new_moduli[0]);
    } else {
        // TODO case
        return RnsPolynomial();
    }
}

} // namespace hehub
