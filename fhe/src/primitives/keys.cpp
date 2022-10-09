#include "keys.h"
#include "common/mod_arith.h"
#include "common/rns_transform.h"

namespace hehub {

RlweKsk::RlweKsk(const RlweSk &sk_new, const RlweSk &sk_orig,
                 const u64 additional_mod) {
    auto sk_new_mult_p_extended = sk_new * additional_mod;
    strict_reduce(sk_new_mult_p_extended); // for more friendly debugging
    // The component corresponding to p after extended is necessarily 0
    sk_new_mult_p_extended.add_components({additional_mod});
    std::fill(sk_new_mult_p_extended.last()->begin(),
              sk_new_mult_p_extended.last()->end(), 0);

    // To encapsulate
    auto extended_moduli = sk_orig.modulus_vec();
    extended_moduli.push_back(additional_mod);
    auto sk_orig_extended(sk_orig);
    intt_negacyclic_inplace_lazy(sk_orig_extended);
    auto extended_part = rns_base_transform(sk_orig_extended, {additional_mod});
    sk_orig_extended.add_components({additional_mod});
    *sk_orig_extended.last() = std::move(extended_part[0]);
    ntt_negacyclic_inplace_lazy(sk_orig_extended);

    auto orig_components = sk_orig.component_count();
    std::vector<std::vector<u64>> rns_composition_basis(orig_components);
    for (size_t i = 0; i < orig_components; i++) {
        rns_composition_basis[i].resize(orig_components + 1, 0);
        rns_composition_basis[i][i] = 1;
    }
    *this = RgswCt(rgsw_encrypt_montgomery(
        sk_orig_extended, sk_new_mult_p_extended, rns_composition_basis));
}

} // namespace hehub
