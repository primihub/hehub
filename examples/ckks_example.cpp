#include "ckks/ckks.h"
#include <cmath>
#include <iostream>

using namespace hehub;

int main() {
    int precision_bits = 30;
    auto params = ckks::create_params(4096, precision_bits);
    CkksSk sk(params);
    auto relin_key = get_relin_key(sk, params.additional_mod);

    auto pt = ckks::encode(1, params);
    auto ct_sum = ckks::encrypt(pt, sk);
    for (int i = 2; i <= 100000; i++) {
        auto pt = ckks::encode(1.0 / i, params);
        auto ct = ckks::encrypt(pt, sk);
        auto ct_squared = ckks::mult(ct, ct, relin_key);
        ct_sum = ckks::add(ct_sum, ct_squared);
    }

    double sum = ckks::decode(ckks::decrypt(ct_sum, sk));
    std::cout << "(" << sum << ", " << M_PI * M_PI / 6 << ")" << std::endl;
}
