#include "fhe/ckks/ckks.h"
#include <cmath>
#include <iostream>
#include <numeric>

using namespace hehub;

int main() {
    int precision_bits = 30;
    auto params = ckks::create_params(4096, precision_bits);
    CkksSk sk(params);
    auto relin_key = get_relin_key(sk, params.additional_mod);

    CkksCt ct_sum;
    for (int i = 1; i <= 10000; i++) {
        auto pt = ckks::encode(1.0 / i, params);
        auto ct = ckks::encrypt(pt, sk);
        auto ct_squared = ckks::mult(ct, ct, relin_key);

        if (i == 1) {
            ct_sum = ct_squared;
        } else {
            ct_sum = ckks::add(ct_sum, ct_squared);
        }
    }

    double sum = ckks::decode(ckks::decrypt(ct_sum, sk));
    std::cout << "(" << sum << ", " << M_PI * M_PI / 6 << ")" << std::endl;
}
