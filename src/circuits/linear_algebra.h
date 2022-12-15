/**
 * @file mat_vec_mul.h
 * @brief Multiplication between matrices and vectors.
 *
 */

#pragma once

#include "ckks/ckks.h"
#include "primitives/keys.h"
#include "range/v3/view/iota.hpp"
#include <vector>

namespace hehub {
namespace ckks {

template <typename T = double,
          typename std::enable_if_t<std::is_same_v<T, double> ||
                                    std::is_same_v<T, cc_double>> * = nullptr>
CkksCt square_matrix_vector_mul(const std::vector<std::vector<T>> &mat,
                                const CkksCt &ct_vec,
                                const std::vector<RotKey> &rot_keys) {
    // a CKKS ciphertext of dimension N has N/2 data slots
    auto vec_dimension = ct_vec[0].dimension() / 2;
    if (mat.size() != vec_dimension) {
        throw std::invalid_argument("Input matrix is not of size (N/2)x(N/2)");
    }
    for (auto &row : mat) {
        if (row.size() != vec_dimension) {
            throw std::invalid_argument(
                "Input matrix is not of size (N/2)x(N/2)");
        }
    }
    for (auto i : ranges::views::ints((size_t)1, vec_dimension)) {
        if (rot_keys[i].step != i) {
            throw std::invalid_argument("Rotation key mismatch");
        }
    }

    auto fast_mod = [=](size_t x) {
        return x & (vec_dimension - 1); // vec_dim ensured to be 2-pow
    };

    CkksCt ct_accumulated;
    CkksParams params = ct_vec[0].params();
    params.initial_scaling_factor = ct_vec.scaling_factor;
    for (auto i : ranges::views::ints((size_t)0, vec_dimension)) {
        std::vector<T> curr_diag(vec_dimension, (T)0);
        for (auto j : ranges::views::ints((size_t)0, vec_dimension)) {
            curr_diag[j] = mat[fast_mod(j + i)][j];
        }
        auto encoded_diag = ckks::simd_encode(curr_diag, params);

        auto ct_prod_diag_vec = ckks::mult_plain(ct_vec, encoded_diag);
        if (i == 0) {
            ct_accumulated = std::move(ct_prod_diag_vec);
        } else {
            auto ct_prod_rotated =
                rotate(ct_prod_diag_vec, rot_keys[i]);
            ct_accumulated = ckks::add(ct_accumulated, ct_prod_rotated);
        }
    }

    return ct_accumulated;
}

} // namespace ckks
} // namespace hehub
