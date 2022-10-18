# HEhub

## Introduction

#### What is homomorphic encryption?

#### Why yet another HE library?
The technology of homomorphic encryption is evolving rapidly, and we aim to provide an easy-to-use, scalable library to facilitate the use and development of homomorphic encryption.

## Building and Installation 
Currently the library only requires header-only third-party dependencies, which require no manual pre-installation. The library is built with CMake (>= 3.14), and tested on Linux with toolchain GCC (>= 7.0) and on MacOS with toolchain Clang (>= 13.0).

To build the library, use the following command to configure and build:
```bash
cmake -S . -B build
cmake --build build
```

To install the library, run the following command:
```bash
sudo cmake --install build
```

## Usage

```c++
#include "ckks/ckks.h"
#include <cmath>
#include <iostream>

using namespace hehub;

int main() {
    int precision_bits = 30;
    CkksParams params(4096, precision_bits);
    CkksSk sk(params);
    auto relin_key = get_relin_key(sk);

    auto pt = ckks::encode(1, params);
    auto ct_sum = ckks::encrypt(pt, sk);
    for (int i = 2; i <= 100000; i++) {
        auto pt = ckks::encode(1.0 / i, params);
        auto ct = ckks::encrypt(pt, sk);
        auto ct_squared = ckks::relinearize(
            ckks::mul(ct, ct), relin_key);
        ct_sum = ckks::add(ct_sum, ct_squared);
    }

    double sum = ckks::decode(ckks::decrypt(ct_sum, sk));
    std::cout << "(" << sum << ", " 
              << M_PI * M_PI / 6 << ")" << std::end;
}

```

## Lisence
