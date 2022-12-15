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

    // prepare the matrix and vector and their product
    const auto VEC_DIM = N / 2;
    MatrixXd mat = MatrixXd::Random(VEC_DIM, VEC_DIM);
    Vector<double, VEC_DIM> vec = Vector<double, VEC_DIM>::Random(VEC_DIM);
    auto prod_vec_ref = mat * vec;

    // convert to STL vector so as to do encoding and encryption
    std::vector<double> vec_copy(vec.data(), vec.data() + vec.size());
    std::vector<std::vector<double>> mat_copy(VEC_DIM);
    for (auto i : ints((size_t)0, VEC_DIM)) {
        MatrixXd row = mat.row(i);
        mat_copy[i] = std::vector<double>(row.data(), row.data() + row.size());
    }

    // encrypt and multiply under encryption
    CkksSk sk(ckks_params);
    std::vector<RotKey> rot_keys(VEC_DIM);
    for (auto step : ints((size_t)1, VEC_DIM)) {
        rot_keys[step] = get_rot_key(sk, ckks_params.additional_mod, step);
    }
    auto vec_ct = ckks::encrypt(ckks::simd_encode(vec_copy, ckks_params), sk);
    auto prod_vec_ct =
        ckks::square_matrix_vector_mul(mat_copy, vec_ct, rot_keys);
    auto prod_vec_decrypted = ckks::simd_decode(ckks::decrypt(prod_vec_ct, sk));

    // check
    const double eps = std::pow(2, -20);
    for (auto i : ints((size_t)0, VEC_DIM)) {
        REQUIRE(std::abs(prod_vec_ref(i) - prod_vec_decrypted[i]) < eps);
    }
}
