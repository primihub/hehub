# HEhub

HEhub is a library for homomorphic encryption and its applications, and is part of the PrimiHub project.

## Introduction

#### What is homomorphic encryption?
Homomorphic encryption (HE) is a cryptographic primitive which allows computation on encrypted messages without the need of decryption or revealing of any infomation about the messages. The notion of fully homomorphic encryption (FHE) is specified by those HE schemes that allow to evaluate circuits of arbituary depth composed with any type of gates, and was first instantiated by Craig Gentry in 2009. Ever since the breakthrough of Gentry, the field of FHE has gone through years of rapid development, and has become a crucial part of state-of-the-art privacy enhancing technologies. In view of the great impact FHE technology has had on cryptography and privacy technology, the authors of papers BV11 and BGV12 about the BGV scheme were awarded Gödel Prize recently.

There are four generations of FHE schemes so far. The schemes of BFV, BGV, and CKKS are currently often applied in leveled mode, in which no bootstrapping procedure will be run and the circuits to evaluate need to be of limited depth. These schemes support SIMD evalutation natively. On the other hand, the scheme of FHEW/TFHE has practical solution of bootstrapping, which enables FHE mode natively. 

#### Why yet another HE library?
There are several open-sourced HE libraries so far. However, the technology of homomorphic encryption is evolving continuously, and we aim to provide an easy-to-use, scalable and efficient library, so as to help developers catch up with the latest development of this field, and to faciliate further research on it. We hope that HEhub can help the cummunity to utilize and explore the future of homomorphic encryption. 

HEhub currently includes the homomorphic encryption schemes BGV, CKKS, and TFHE, etc., and will further feature various schemes, frequent circuits and application interface on homomorphic encryption. As part of the PrimiHub project, HEhub is an essential tool helping us explore the field of privacy enhancing technologies.

## Building and Installation 
Currently the library only requires header-only third-party dependencies, which need no manual pre-installation. The library is built with CMake (>= 3.14), and tested on Linux with toolchain GCC (>= 7.0) and on MacOS with toolchain Clang (>= 13.0).

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

The usage of homomorphic encryption with HEhub is very simple. Below is an example computing the Basel series with the CKKS scheme.

```cpp
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
        auto ct_squared = ckks::mult(ct, ct, relin_key);
        ct_sum = ckks::add(ct_sum, ct_squared);
    }

    double sum = ckks::decode(ckks::decrypt(ct_sum, sk));
    std::cout << "(" << sum << ", " 
              << M_PI * M_PI / 6 << ")" << std::end;
}

```

## Contributing
We strongly welcome developers to contribute code for HEhub. On how to participate, see [PrimiHub Open Source Community Governance](http://docs.primihub.com/docs/primihub-community).
