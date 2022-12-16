#include "Eigen/Dense"
#include "catch2/catch.hpp"
#include "linear_algebra.h"
#include "range/v3/view/iota.hpp"

using namespace hehub;
using namespace Eigen;
using namespace ranges::views;

TEST_CASE("matrix vector mul", "[.]") {
    // CKKS parameters
    const size_t N = 4096;
    auto ckks_params = ckks::create_params(N, 35);
    SECTION("square matrix") {
        // let the vector length be right equal to slot count
        // and let the matrix be square
        const auto vec_dim = N / 2;
        MatrixXd mat = MatrixXd::Random(vec_dim, vec_dim);
        Vector<double, vec_dim> vec = Vector<double, vec_dim>::Random(vec_dim);
        auto prod_vec_ref = mat * vec;

        // convert to STL vector so as to do encoding and encryption
        std::vector<double> vec_copy(vec.data(), vec.data() + vec.size());
        std::vector<std::vector<double>> mat_copy(vec_dim);
        for (auto i : ints((size_t)0, vec_dim)) {
            MatrixXd row = mat.row(i);
            mat_copy[i] =
                std::vector<double>(row.data(), row.data() + row.size());
        }

        // encrypt and multiply under encryption
        CkksSk sk(ckks_params);
        auto slot_count = N / 2;
        std::vector<RotKey> rot_keys(slot_count);
        // only need rotation key for 1 step
        rot_keys[1] = get_rot_key(sk, ckks_params.additional_mod, 1);
        auto vec_ct =
            ckks::encrypt(ckks::simd_encode(vec_copy, ckks_params), sk);
        auto prod_vec_ct =
            ckks::matrix_vector_mul_short(mat_copy, vec_ct, rot_keys);
        auto prod_vec_decrypted =
            ckks::simd_decode(ckks::decrypt(prod_vec_ct, sk));

        // check
        const double eps = std::pow(2, -20);
        for (auto i : ints((size_t)0, vec_dim)) {
            REQUIRE(std::abs(prod_vec_ref(i) - prod_vec_decrypted[i]) < eps);
        }
    }
    SECTION("small matrix") {
        // let the vector length be smaller than slot count
        // and let the matrix be short in shape
        const auto vec_dim = N / 4;
        const auto mat_height = vec_dim - 1;
        const auto mat_width = vec_dim;
        MatrixXd mat = MatrixXd::Random(mat_height, mat_width);
        Vector<double, vec_dim> vec = Vector<double, vec_dim>::Random(vec_dim);
        auto prod_vec_ref = mat * vec;

        // convert to STL vector so as to do encoding and encryption
        std::vector<double> vec_copy(vec.data(), vec.data() + vec.size());
        std::vector<std::vector<double>> mat_copy(mat_height);
        for (auto i : ints((size_t)0, mat_height)) {
            MatrixXd row = mat.row(i);
            mat_copy[i] =
                std::vector<double>(row.data(), row.data() + row.size());
        }

        // encrypt and multiply under encryption
        CkksSk sk(ckks_params);
        auto slot_count = N / 2;
        std::vector<RotKey> rot_keys(slot_count);
        // need (2*mat_width)-many rotation keys
        const auto steps = ckks::mv_mul_requiring_steps(slot_count, mat_width);
        for (auto step : steps) {
            rot_keys[step] = get_rot_key(sk, ckks_params.additional_mod, step);
        }
        auto vec_ct =
            ckks::encrypt(ckks::simd_encode(vec_copy, ckks_params), sk);
        auto prod_vec_ct =
            ckks::matrix_vector_mul_short(mat_copy, vec_ct, rot_keys);
        auto prod_vec_decrypted =
            ckks::simd_decode(ckks::decrypt(prod_vec_ct, sk));

        // check
        const double eps = std::pow(2, -20);
        for (auto i : ints((size_t)0, mat_height)) {
            REQUIRE(std::abs(prod_vec_ref(i) - prod_vec_decrypted[i]) < eps);
        }
    }
}
