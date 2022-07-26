/**
 * @file type_defs.h
 * @brief Define short alias for common integer types.
 *
 */
#pragma once
#include <stdint.h>

namespace hehub {

using i8 = __int8_t;
using i16 = __int16_t;
using i32 = __int32_t;
using i64 = __int64_t;
using i128 = __int128_t;

using u8 = __uint8_t;
using u16 = __uint16_t;
using u32 = __uint32_t;
using u64 = __uint64_t;
using u128 = __uint128_t;

} // namespace hehub
