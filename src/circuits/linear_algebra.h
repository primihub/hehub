/**
 * @file linear_algebra.h
 * @brief Linear algebra functionalities.
 *
 */

#pragma once

#include "fhe/ckks/ckks.h"
#include "fhe/common/ntt.h"
#include "fhe/primitives/keys.h"
#include "range/v3/view/iota.hpp"
#include <vector>

namespace hehub {
namespace ckks {

/**
 * @brief Helper function which tells which rotation steps are needed for
 * matrix-vector multiplication. (Only used in the case when matrix width < slot
 * count.)
 * @param slot_count The number of slots in a CKKS ciphertext.
 * @param matrix_width Width of the matrix (i.e., the dimension of the vector to
 * be multiplied).
 * @return std::vector<size_t>
 */
std::vector<size_t> mv_mul_requiring_steps(size_t slot_count,
                                           size_t matrix_width) {
    std::vector<size_t> steps;
    for (auto step : ranges::views::ints((size_t)1, matrix_width)) {
        steps.push_back(step);
        steps.push_back(step + slot_count - matrix_width);
    }
    return steps;
}

/**
 * @brief Multiplication between an encrypted vector and a matrix in
 * cleartext, outputting the product vector under encryption. The input
 * vector is encrypted as a single CKKS ciphertext, and the matrix is short
 * (i.e. matrix height <= slot count).
 * @tparam T Double or complex double
 * @param mat The input matrix
 * @param ct_vec The ciphertext carrying the input vector
 * @param rot_key_set Rotation key set (currently requires keys for steps 1
 * ~ matrix width)
 * @return CkksCt
 */
template <typename T = double,
          typename std::enable_if_t<std::is_same_v<T, double> ||
                                    std::is_same_v<T, cc_double>> * = nullptr>
CkksCt matrix_vector_mul_short(const std::vector<std::vector<T>> &mat,
                               const CkksCt &ct_vec,
                               const std::vector<RotKey> &rot_key_set) {
    // a CKKS ciphertext of dimension N has N/2 data slots
    auto slot_count = ct_vec[0].dimension() / 2;
    auto matrix_height = mat.size();
    if (matrix_height > slot_count) {
        throw std::invalid_argument(
            "Input matrix is too tall.\nHint: invoke the "
            "version for tall matrix.");
    }
    size_t matrix_width = 0;
    for (auto &row : mat) {
        matrix_width = std::max(row.size(), matrix_width);
    }
    if (matrix_width > slot_count) {
        throw std::invalid_argument(
            "Input matrix is too wide. Choose larger CKKS parameters.");
    }
    if (matrix_width > slot_count / 2 && matrix_width != slot_count) {
        throw std::invalid_argument(
            "Input matrix's width should be equal to slot count, or not "
            "greater than half of slot count. Choose proper CKKS parameters.");
    }
    if (rot_key_set.size() != slot_count) {
        throw std::invalid_argument("Rotation key set is malformed.");
    }

    const bool full_width = (matrix_width == slot_count);
    const auto steps = full_width
                           ? std::vector<size_t>{1}
                           : mv_mul_requiring_steps(slot_count, matrix_width);
    if (full_width) {
        if (rot_key_set[1].empty()) {
            throw std::invalid_argument(
                "Required rotation key not generated. (Need key for 1 "
                "step.)");
        }
    } else {
        for (auto step : steps) {
            if (rot_key_set[step].empty()) {
                throw std::invalid_argument(
                    "Required rotation key not generated. (Need key for " +
                    std::to_string(step) + " steps.)");
            }
        }
    }

    auto ct_vec_rotating(ct_vec);
    CkksCt ct_accumulated;
    CkksParams params = ct_vec[0].params();
    params.initial_scaling_factor = ct_vec.scaling_factor;
    for (auto i : ranges::views::ints((size_t)0, matrix_width)) {
        std::vector<T> curr_diag(slot_count, (T)0);
        for (auto j : ranges::views::ints((size_t)0, matrix_height)) {
            curr_diag[j] = mat[j][(j + matrix_width - i) % matrix_width];
        }
        auto encoded_diag = simd_encode(curr_diag, params);

        auto ct_prod_diag_vec = mult_plain(ct_vec_rotating, encoded_diag);

        if (i == 0) {
            ct_accumulated = std::move(ct_prod_diag_vec);
        } else {
            ct_accumulated = add(ct_accumulated, ct_prod_diag_vec);
        }

        if (i != matrix_width - 1) {
            // update the encryted vector by rotation
            if (full_width) {
                ct_vec_rotating = rotate(ct_vec_rotating, rot_key_set[1]);
            } else {
                auto next_step = i + 1;
                ct_vec_rotating = add(
                    rotate(ct_vec, rot_key_set[next_step]),
                    rotate(ct_vec,
                           rot_key_set[next_step + slot_count - matrix_width]));
            }
        }
    }

    rescale_inplace(ct_accumulated);
    return ct_accumulated;
}

} // namespace ckks
} // namespace hehub
