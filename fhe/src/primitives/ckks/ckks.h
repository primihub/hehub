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
 * 
 */
struct CkksQuadraticCt : public std::array<RnsPolynomial, 3> {
    using std::array<RnsPolynomial, 3>::array;

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
    static std::vector<T> simd_decode(const CkksPt &pt, size_t data_size = 0);

    /**
     * @brief TODO
     *
     * @param pt
     * @param sk
     * @return CkksCt
     */
    inline static CkksCt encrypt(const CkksPt &pt, const RlweSk &sk) {
        CkksCt ct = ::hehub::encrypt_core(pt, sk);
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
    inline static CkksPt decrypt(const CkksCt &ct, const RlweSk &sk) {
        CkksPt pt = ::hehub::decrypt_core(ct, sk);
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
    static CkksCt add(const CkksCt &ct1, const CkksCt &ct2);

    /**
     * @brief TODO
     * 
     * @param ct 
     * @param pt 
     * @return CkksCt 
     */
    static CkksCt add_plain(const CkksCt &ct, const CkksPt &pt);

    /**
     * @brief TODO
     * 
     * @param ct1 
     * @param ct2 
     * @return CkksCt 
     */
    static CkksCt sub(const CkksCt &ct1, const CkksCt &ct2);

    /**
     * @brief TODO
     * 
     * @param ct 
     * @param pt 
     * @return CkksCt 
     */
    static CkksCt sub_plain(const CkksCt &ct, const CkksPt &pt);

    /**
     * @brief TODO
     * 
     * @param ct 
     * @param pt 
     * @return CkksCt 
     */
    static CkksCt mult_plain(const CkksCt &ct, const CkksPt &pt);

    /**
     * @brief TODO
     * 
     * @param ct1 
     * @param ct2 
     * @return CkksQuadraticCt 
     */
    static CkksQuadraticCt mult_low_level(const CkksCt &ct1, const CkksCt &ct2);

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
