#pragma once
#include "type_defs.h"

namespace hehub
{

inline u64 __bit_rev_naive(u64 x, const size_t bit_len) {
    if (bit_len <= 1) { return x; }
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

} // namespace hehub
