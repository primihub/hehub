/**
 * @file ckks.h
 * @brief TODO
 *
 */

#pragma once

#include "common/type_defs.h"
#include "primitives/rlwe.h"
#include <complex>

namespace hehub {

/**
 * @brief TODO
 * 
 */
struct CkksPt: public RlwePt {
    using RlwePt::RlwePt;

    /// @brief TODO
    double scaling_factor = 1.0;
};

/**
 * @brief TODO
 * 
 */
struct CkksCt: public RlweCt {
    using RlweCt::RlweCt;

    /// @brief TODO
    double scaling_factor = 1.0;
};

struct ckks {
    /**
     * @brief TODO
     *
     * @param data
     * @param scaling_factor
     * @param pt_poly_dim
     * @return CkksPt
     */
    static CkksPt simd_encode_cc(const std::vector<cc_double> &data,
                                 const double scaling_factor,
                                 const PolyDimensions &pt_poly_dim);

    /**
     * @brief TODO
     *
     * @param data
     * @param scaling_factor
     * @param pt_poly_dim
     * @return CkksPt
     */
    static CkksPt simd_encode(const std::vector<cc_double> &data,
                              const double scaling_factor,
                              const PolyDimensions &pt_poly_dim) {
        return simd_encode_cc(data, scaling_factor, pt_poly_dim);
    }

    /**
     * @brief TODO
     *
     * @param data
     * @param scaling_factor
     * @param pt_poly_dim
     * @return CkksPt
     */
    static CkksPt simd_encode(const std::vector<double> &data,
                              const double scaling_factor,
                              const PolyDimensions &pt_poly_dim) {
        std::vector<cc_double> data_cc;
        for (auto d : data) {
            data_cc.push_back(cc_double(d));
        }
        return simd_encode_cc(data_cc, scaling_factor, pt_poly_dim);
    }

    /**
     * @brief TODO
     *
     * @param pt
     * @param data_size
     * @return std::vector<cc_double>
     */
    static std::vector<cc_double> simd_decode_cc(const CkksPt &pt,
                                                 size_t data_size = 0);

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
                                      std::is_same<T, cc_double>::value>::type
                  * = nullptr>
    static std::vector<T> simd_decode(const CkksPt &pt,
                                      size_t data_size = 0);

    /**
     * @brief TODO
     * 
     * @param ct 
     * @param dropping_primes 
     */
    static void rescale_inplace(CkksCt &ct, size_t dropping_primes = 1);

private:
    // Instantiation is diabled.
    ckks();
};

} // namespace hehub
