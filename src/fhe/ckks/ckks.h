/**
 * @file ckks.h
 * @brief TODO
 *
 */

#pragma once

#include "common/type_defs.h"
#include "primitives/keys.h"
#include "primitives/rgsw.h"
#include "primitives/rlwe.h"
#include <complex>

namespace hehub {
namespace ckks {

struct CkksParams : public RlweParams {
    CkksParams(RlweParams &&other) : RlweParams(std::move(other)) {}

    u64 additional_mod = 1;

    double initial_scaling_factor = 1.0;
};

/**
 * @brief Create a CkkaParams object, which determines the dimension, moduli and
 * initial scaling factor from input arguments.
 * @param dimension The RLWE dimension.
 * @param moduli_bits A list of integers specifying bit length of each modulus
 * of the RNS.
 * @param additional_mod_bits An integer specifying bit length of the additional
 * modulus.
 * @param initial_scaling_factor The initial scaling factor.
 * @return CkksParams
 */
CkksParams create_params(size_t dimension, std::vector<int> moduli_bits,
                         int additional_mod_bits,
                         double initial_scaling_factor);

/**
 * @brief Create a params object, which automatically determines the dimension,
 * moduli and initial scaling factor.
 * @param dimension The RLWE dimension.
 * @param initial_scaling_factor The initial scaling factor.
 * @return CkksParams
 */
CkksParams create_params(size_t dimension, double initial_scaling_factor);

/**
 * @brief TODO
 *
 */
struct CkksPt : public RlwePt {
    using RlwePt::RlwePt;

    /// @brief TODO
    /// @param other
    CkksPt(RlwePt &&other) : RlwePt(std::move(other)) {}

    /// @brief TODO
    double scaling_factor = 1.0;
};

/**
 * @brief TODO
 *
 */
struct CkksCt : public RlweCt {
    using RlweCt::RlweCt;

    /// @brief TODO
    /// @param other
    CkksCt(RlweCt &&other) : RlweCt(std::move(other)) {}

    /// @brief TODO
    double scaling_factor = 1.0;
};

/**
 * @brief TODO
 * (TODO: merge this into CkksCt)
 */
struct CkksQuadraticCt : public std::array<RnsPolynomial, 3> {
    using std::array<RnsPolynomial, 3>::array;

    /// @brief TODO
    double scaling_factor = 1.0;
};

/**
 * @brief TODO
 *
 * @param data
 * @param scaling_factor
 * @param pt_params
 * @return CkksPt
 */
CkksPt simd_encode(const std::vector<cc_double> &data,
                   const double scaling_factor, const CkksParams &pt_params);

/**
 * @brief TODO
 *
 * @param data
 * @param scaling_factor
 * @param pt_params
 * @return CkksPt
 */
CkksPt simd_encode(const std::vector<double> &data, const double scaling_factor,
                   const CkksParams &pt_params);

/**
 * @brief TODO
 *
 * @param datum
 * @param scaling_factor
 * @param pt_params
 * @return CkksPt
 */
inline CkksPt encode(const cc_double datum, const double scaling_factor,
                     const CkksParams &pt_params) {
    std::vector datum_rep(pt_params.dimension / 2, datum);
    return simd_encode(datum_rep, scaling_factor, pt_params);
}

/**
 * @brief TODO
 *
 * @param datum
 * @param scaling_factor
 * @param pt_params
 * @return CkksPt
 */
inline CkksPt encode(const double datum, const double scaling_factor,
                     const CkksParams &pt_params) {
    std::vector datum_rep(pt_params.dimension / 2, datum);
    return simd_encode(datum_rep, scaling_factor, pt_params);
}

/**
 * @brief TODO
 *
 * @tparam T
 * @param pt
 * @param data_size
 * @return std::vector<T>
 */
template <typename T = double,
          typename std::enable_if<std::is_same<T, double>::value ||
                                  std::is_same<T, cc_double>::value>::type * =
              nullptr>
std::vector<T> simd_decode(const CkksPt &pt, size_t data_size = 0);

/**
 * @brief TODO
 *
 * @param pt
 * @param sk
 * @return CkksCt
 */
inline CkksCt encrypt(const CkksPt &pt, const RlweSk &sk) {
    CkksCt ct = encrypt_core(pt, sk);
    ct.scaling_factor = pt.scaling_factor;
    return ct;
}

/**
 * @brief TODO
 *
 * @param ct
 * @param sk
 * @return CkksCt
 */
inline CkksPt decrypt(const CkksCt &ct, const RlweSk &sk) {
    CkksPt pt = decrypt_core(ct, sk);
    pt.scaling_factor = ct.scaling_factor;
    return pt;
}

/**
 * @brief TODO
 *
 * @param ct1
 * @param ct2
 * @return CkksCt
 */
CkksCt add(const CkksCt &ct1, const CkksCt &ct2);

/**
 * @brief TODO
 *
 * @param ct
 * @param pt
 * @return CkksCt
 */
CkksCt add_plain(const CkksCt &ct, const CkksPt &pt);

/**
 * @brief TODO
 *
 * @param ct1
 * @param ct2
 * @return CkksCt
 */
CkksCt sub(const CkksCt &ct1, const CkksCt &ct2);

/**
 * @brief TODO
 *
 * @param ct
 * @param pt
 * @return CkksCt
 */
CkksCt sub_plain(const CkksCt &ct, const CkksPt &pt);

/**
 * @brief TODO
 *
 * @param ct
 * @param pt
 * @return CkksCt
 */
CkksCt mult_plain(const CkksCt &ct, const CkksPt &pt);

/**
 * @brief TODO
 *
 * @param ct1
 * @param ct2
 * @return CkksQuadraticCt
 */
CkksQuadraticCt mult_low_level(const CkksCt &ct1, const CkksCt &ct2);

/**
 * @brief TODO
 *
 * @param ct
 * @param relin_key
 * @return CkksCt
 */
CkksCt relinearize(const CkksQuadraticCt &ct, const RlweKsk &relin_key);

/**
 * @brief TODO
 *
 * @param ct1
 * @param ct2
 * @param relin_key
 * @return CkksCt
 */
inline CkksCt mult(const CkksCt &ct1, const CkksCt &ct2,
                   const RlweKsk &relin_key) {
    auto ct_prod = mult_low_level(ct1, ct2);
    return relinearize(ct_prod, relin_key);
}

/**
 * @brief TODO
 *
 * @param ct
 * @param conj_key
 * @return CkksCt
 */
CkksCt conjugate(const CkksCt &ct, const RlweKsk &conj_key);

/**
 * @brief TODO
 *
 * @param ct
 * @param rot_key
 * @param step
 * @return CkksCt
 */
CkksCt rotate(const CkksCt &ct, const RlweKsk &rot_key, const size_t step);

/**
 * @brief TODO
 *
 * @param ct
 * @param dropping_primes
 */
void rescale_inplace(CkksCt &ct, size_t dropping_primes = 1);

} // namespace ckks
} // namespace hehub

// export types
namespace hehub {

using CkksParams = ckks::CkksParams;

using CkksSk = RlweSk;

using CkksPt = ckks::CkksPt;

using CkksCt = ckks::CkksCt;

} // namespace hehub
