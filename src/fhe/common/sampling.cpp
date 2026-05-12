#include "sampling.h"
#include "ntt.h"
#include "range/v3/view/zip.hpp"
#include "rns.h"
#include "type_defs.h"
#include <cmath>
#include <random>

using namespace ranges::views;

namespace hehub {

// ── Fast thread-local xorshift1024* PRNG ──────────────────────────────
// Not meant for cryptographic use; sufficient for HE noise distributions.
struct alignas(64) XorShift1024 {
    u64 s[16];
    int p;

    XorShift1024() {
        std::random_device rd;
        for (int i = 0; i < 16; i++) {
            s[i] = ((u64)rd() << 32) | rd();
            if (s[i] == 0) s[i] = 1;
        }
        p = 0;
    }

    u64 next() {
        u64 s0 = s[p];
        u64 s1 = s[p = (p + 1) & 15];
        s1 ^= s1 << 31;
        s[p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30);
        return s[p] * 1181783497276652981ULL;
    }

    u64 next_bounded(u64 limit) {
        u64 t = -limit % limit;
        u64 r;
        while ((r = next()) < t);
        return r % limit;
    }

    double next_double() {
        return (next() >> 11) * (1.0 / 9007199254740992.0);
    }
};

static thread_local XorShift1024 tl_rng;

RnsPolynomial get_rand_ternary_poly(const RnsPolyParams &params) {
    RnsPolynomial tern_poly(params);
    auto dimension = params.dimension;

    std::vector<i8> ternary_integers(dimension);
    for (auto &t : ternary_integers) {
        // Generate {-1, 0, 1} uniformly
        u64 r = tl_rng.next_bounded(3);
        t = (i8)(r - 1);
    }

    for (auto [component, modulus] : zip(tern_poly, tern_poly.modulus_vec())) {
        for (auto [coeff, ternary_int] : zip(component, ternary_integers)) {
            coeff = modulus + (u64)ternary_int;
            coeff -= (coeff >= modulus) ? modulus : 0;
        }
    }

    ntt_negacyclic_inplace_lazy(tern_poly);
    return tern_poly;
}

RnsPolynomial get_rand_uniform_poly(const RnsPolyParams &params,
                                    PolyRepForm form) {
    auto dimension = params.dimension;
    RnsPolynomial rand_rns_poly(params);

    for (auto [component, modulus] :
         zip(rand_rns_poly, rand_rns_poly.modulus_vec())) {
        for (auto &coeff : component) {
            coeff = tl_rng.next_bounded(modulus);
        }
    }

    rand_rns_poly.rep_form = form;
    return rand_rns_poly;
}

/// Box-Muller transform: generates two independent standard normal samples.
static void box_muller_std_norm(double &z1, double &z2) {
    double u1, u2, r;
    do {
        u1 = tl_rng.next_double() * 2.0 - 1.0;
        u2 = tl_rng.next_double() * 2.0 - 1.0;
        r = u1 * u1 + u2 * u2;
    } while (r >= 1.0 || r == 0.0);

    double c = std::sqrt(-2.0 * std::log(r) / r);
    z1 = u1 * c;
    z2 = u2 * c;
}

RnsPolynomial get_rand_gaussian_poly(const RnsPolyParams &params,
                                     double std_dev) {
    auto dimension = params.dimension;
    RnsPolynomial gaussian_poly(params);

    auto bound = std_dev * 6;

    std::vector<double> gaussians(dimension);
    for (size_t i = 0; i < dimension; i += 2) {
        double z1, z2;
        do {
            box_muller_std_norm(z1, z2);
            z1 *= std_dev;
            z2 *= std_dev;
        } while (std::abs(z1) > bound || (i + 1 < dimension && std::abs(z2) > bound));

        gaussians[i] = z1;
        if (i + 1 < dimension) {
            gaussians[i + 1] = z2;
        }
    }

    for (auto [component, modulus] :
         zip(gaussian_poly, gaussian_poly.modulus_vec())) {
        for (auto [coeff, gaussian] : zip(component, gaussians)) {
            coeff = modulus + (u64)std::llround(gaussian);
            coeff -= (coeff >= modulus) ? modulus : 0;
        }
    }

    ntt_negacyclic_inplace_lazy(gaussian_poly);
    return gaussian_poly;
}

RnsPolynomial get_zero_poly(const RnsPolyParams &params, PolyRepForm form) {
    RnsPolynomial rns_poly(params);
    rns_poly.rep_form = form;
    for (auto &component : rns_poly) {
        std::fill(component.begin(), component.end(), 0);
    }
    return rns_poly;
}

} // namespace hehub
