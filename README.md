# HEhub

HEhub is a library for homomorphic encryption and its applications, and is part of the [PrimiHub](https://github.com/primihub/primihub) project.

## Introduction

#### What is homomorphic encryption?
Homomorphic encryption (HE) is a cryptographic primitive which allows computation on encrypted messages without the need of decryption or revealing of any infomation about the messages. The notion of fully homomorphic encryption (FHE) is specified by those HE schemes that allow to evaluate circuits of arbituary depth composed with any type of gates, and was first instantiated by Craig Gentry in 2009. Ever since the breakthrough of Gentry, the field of FHE has gone through years of rapid development, and has become a crucial part of state-of-the-art privacy enhancing technologies. In view of the great impact FHE technology has had on cryptography and privacy technology, the authors of papers BV11 and BGV12 about the BGV scheme were awarded Gödel Prize recently.

There are four generations of FHE schemes so far. The schemes of BFV, BGV, and CKKS are currently often applied in leveled mode, in which no bootstrapping procedure will be run and the circuits to evaluate need to be of limited depth. These schemes support SIMD evalutation natively. On the other hand, the scheme of FHEW/TFHE has practical solution of bootstrapping, which enables FHE mode natively. 

#### Why yet another HE library?
There are several open-sourced HE libraries so far. However, the technology of homomorphic encryption is evolving continuously, and we aim to provide an easy-to-use, scalable and efficient library, so as to help developers catch up with the latest development of this field, and to faciliate further research on it. We hope that HEhub can help the cummunity to utilize and explore the future of homomorphic encryption. 

HEhub currently includes the homomorphic encryption schemes BGV, CKKS, and TFHE, etc., and will further feature various schemes, frequent circuits and application interface on homomorphic encryption. As part of the PrimiHub project, HEhub is an essential tool helping us explore the field of privacy enhancing technologies.

## Building and Installation 
Currently the library only requires header-only third-party dependencies, which need no manual pre-installation. The library is built with CMake (>= 3.14), and tested on Linux with toolchain GCC (>= 8.0) and on MacOS with toolchain Clang (>= 12.0).

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
    auto params = ckks::create_params(4096, precision_bits);
    CkksSk sk(params);
    auto relin_key = get_relin_key(sk, params.additional_mod);

    CkksCt ct_sum;
    for (int i = 1; i <= 100000; i++) {
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

```

## Benchmarks
We tested the performance of HEhub compiled with Clang-12.0.5 and run on an Intel i7-9750H @ 2.60GHz. _Note: The code for benchmark is still incomplete since our effort is limited currently. We will list more benchmark results later._

| parameter set |  NTT  |  INTT  | CKKS<br>encode +<br>encrypt | CKKS<br>decrypt +<br>decode |
| ------------- |  ---  |  ----  | --------------------------- | --------------------------- |
| N = 1024      |  7 us |   9 us |                             |                             |
| N = 2048      | 14 us |  19 us |                             |                             |
| N = 4096      | 30 us |  37 us |                      426 us |                      237 us |
| N = 8192      | 68 us |  85 us |                    1.730 ms |                      842 us |
| N = 16384     | 142 us| 195 us |                    6.776 ms |                    3.824 ms |
| N = 32768     | 330 us| 406 us |                   27.414 ms |                   18.623 ms |

## How to contribute
If you want to contribute to this project, feel free to create an issue at our [Issue](https://github.com/primihub/primihub/issues) page (e.g., documentation, new idea and proposal).

Also, you can learn about our community [PrimiHub Open Source Community Governance](http://docs.primihub.com/docs/primihub-community)

This is an active open source project for everyone, and we are always open to everyone who want to use this system or contribute to it.

## Community
* Slack: [primihub.slack.com](https://join.slack.com/t/primihub/shared_invite/zt-1iftyi7x0-n_HqllTgPfoEcgqw5UzoYw)
* Wechat Official Account:

![wechat_helper](./doc/wechat.jpeg)
