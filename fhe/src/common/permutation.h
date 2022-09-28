#pragma once
#include "type_defs.h"

namespace hehub {

inline u64 __bit_rev_naive(u64 x, int bit_len) {
#ifdef HEHUB_DEBUG
    if (bit_len < 0 || bit_len > 64) {
        throw std::invalid_argument("bit_len");
    }
    if (bit_len != 64 && x >= (1ULL << bit_len)) {
        throw std::invalid_argument("x");
    }
#endif
    if (bit_len <= 1) {
        return x;
    }
    u64 mh = (1 << (bit_len - 1));
    u64 ml = 1;
    for (; mh >= ml; mh >>= 1, ml <<= 1) {
        if ((mh & x) && (ml & x)) {
        } else if (!(mh & x) && !(ml & x)) {
        } else if (!(mh & x) && (ml & x)) {
            x |= mh;
            x ^= ml;
        } else {
            x |= ml;
            x ^= mh;
        }
    }
    return x;
}

inline u64 __bit_rev_naive_16(u64 x, int bit_len) {
#ifdef HEHUB_DEBUG
    if (bit_len < 0 || bit_len > 16) {
        throw std::invalid_argument("bit_len");
    }
    if (x >= (1ULL << bit_len)) {
        throw std::invalid_argument("x");
    }
#endif
    x = ((x & 0xFF00FF00) >> 8) | ((x & 0x00FF00FF) << 8);
    x = ((x & 0xF0F0F0F0) >> 4) | ((x & 0x0F0F0F0F) << 4);
    x = ((x & 0xCCCCCCCC) >> 2) | ((x & 0x33333333) << 2);
    x = ((x & 0xAAAAAAAA) >> 1) | ((x & 0x55555555) << 1);
    return x >> (16 - bit_len);
}

} // namespace hehub
