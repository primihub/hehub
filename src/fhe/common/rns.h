/**
 * @file rns.h
 * @brief TODO
 *
 */
#pragma once

#include "type_defs.h"
#include <sstream>
#include <vector>

namespace hehub {

/**
 * @class RnsIntVec
 * @brief An RnsIntVec represents an element of Z_q[X]/(X^n + 1), where q is
 * a big composite integer, and n is a certain power of 2. The object contains a
 * group of polynomials whose coefficients of terms of same degree constitute an
 * residue number system (RNS) whose basis is the co-prime factors of q. The
 * group of polynomials are each represented by n coefficients (or n values from
 * the transformation of NTT). The X^n + 1 is sometimes called the
 * polynomial modulus of this object.
 */
class RnsIntVec {
public:
    /**
     * @brief Specifying the dimensions of an RnsIntVec, i.e. the length of
     * polynomial, the number of RNS components and the set of moduli. Used for
     * initializing an RnsIntVec object.
     */
    struct Params {
        /// The length of polynomial.
        size_t dimension = 0;

        /// The number of RNS components.
        size_t component_count;

        /// The set of RNS moduli.
        std::vector<u64> moduli;
    };

    /**
     * @class ComponentData
     * @brief A dynamic array for storing all coefficients of a component
     * polynomial of an RnsIntVec object.
     */
    class ComponentData {
    public:
        /**
         * @brief Construct a new Simple Poly object.
         */
        ComponentData() {}

        /**
         * @brief Construct a new ComponentData object.
         * @param dimension The length of the polynomial.
         */
        ComponentData(const size_t dimension);

        /**
         * @brief Create a ComponentData object by copying from another one.
         * @param other The source ComponentData object.
         */
        ComponentData(const ComponentData &other);

        /**
         * @brief Move another ComponentData object to construct this.
         * @param other The source ComponentData object.
         */
        ComponentData(ComponentData &&other) noexcept {
            *this = std::move(other);
        }

        /**
         * @brief Destroy the ComponentData object.
         */
        ~ComponentData();

        /**
         * @brief Returns the i-th coefficient.
         * @param idx The index.
         * @return u64
         */
        inline u64 &operator[](const int idx) { return data_[idx]; }

        /**
         * @brief Returns the i-th coefficient.
         * @param idx The index.
         * @return const u64
         */
        inline const u64 operator[](const int idx) const { return data_[idx]; }

        /**
         * @brief Copies a given object to the current one.
         * @param copying The given Simple Poly object.
         * @return ComponentData &
         */
        ComponentData &operator=(const ComponentData &copying);

        /**
         * @brief Moves a given object to the current one.
         * @param moving The given Simple Poly object.
         * @return ComponentData &
         */
        ComponentData &operator=(ComponentData &&moving) noexcept;

        /**
         * @brief Returns whether the input object is equivalent to this object.
         * @param other The input comparing object.
         * @return bool
         */
        inline bool operator==(const ComponentData &other) const {
            if (dimension_ != other.dimension_)
                return false;
            for (int i = 0; i < dimension_; i++) {
                if ((*this)[i] != other[i]) {
                    return false;
                }
            }
            return true;
        }

        /**
         * @brief Returns whether the input object is different from this
         * object.
         * @param comparing The input comparing object.
         * @return bool
         */
        inline bool operator!=(const ComponentData &comparing) const {
            return !((*this) == comparing);
        }

        /**
         * @brief Returns the data as a non-const pointer.
         * @return u64 *
         */
        inline u64 *data() { return data_; }

        /**
         * @brief Returns the data as a const pointer.
         * @return const u64 *
         */
        inline const u64 *data() const { return data_; }

        /**
         * @brief Returns the begin of data as a non-const pointer.
         * @return u64 *
         */
        inline u64 *begin() { return data_; }

        /**
         * @brief Returns the begin of data as a const pointer.
         * @return const u64 *
         */
        inline const u64 *begin() const { return data_; }

        /**
         * @brief Returns the end of data as a non-const pointer.
         * @return u64 *
         */
        inline u64 *end() { return data_ + dimension_; }

        /**
         * @brief Returns the end of data as a const pointer.
         * @return const u64 *
         */
        inline const u64 *end() const { return data_ + dimension_; }

    private:
        /// Pointer as a dynamic array for storing the coefficients.
        u64 *data_ = nullptr;

        /// The length of the polynomial, i.e. the number of coefficients or NTT
        /// values.
        size_t dimension_ = 0;
    };

    /**
     * @brief Representation form of the polynomial.
     */
    enum class RepForm { coeff, value };

    /**
     * @brief Construct an empty RnsIntVec object.
     */
    RnsIntVec() {}

    /**
     * @brief Construct a new RnsIntVec object.
     * @param dimension The length of each component polynomial, i.e. the number
     * of coefficients (or NTT values). This should be a power of 2.
     * @param components The number of components, i.e. polynomials each
     * modulo different integer.
     * @param moduli The modulus set of the RNS of coefficients.
     */
    RnsIntVec(const size_t dimension, const size_t components,
              const std::vector<u64> &moduli);

    /**
     * @brief Construct a new RnsIntVec object.
     * @param params A parameter set specifying the length of polynomial, the
     * number of RNS components and the modulus set.
     */
    RnsIntVec(const Params &params);

    /**
     * @brief Creates an RnsIntVec object by copying from another one.
     * @param other The source RnsIntVec object.
     */
    RnsIntVec(const RnsIntVec &other) = default;

    /**
     * @brief Creates an RnsIntVec object by moving another one.
     * @param other The source RnsIntVec object.
     */
    RnsIntVec(RnsIntVec &&other) = default;

    /**
     * @brief Copies a given RnsIntVec to this.
     * @param other The given RnsIntVec object.
     * @return RnsIntVec &
     */
    RnsIntVec &operator=(const RnsIntVec &other) = default;

    /**
     * @brief Moves a given RnsIntVec to replace this.
     * @param other The given RnsIntVec object.
     * @return RnsIntVec &
     */
    RnsIntVec &operator=(RnsIntVec &&other) = default;

    /**
     * @brief Compares this with another RnsIntVec, returning true only when
     * all spec and data coincide.
     * @param other The RnsIntVec object for comparison.
     * @return const bool
     */
    inline const bool operator==(const RnsIntVec &other) const {
        return log_dimension_ == other.log_dimension_ &&
               dimension_ == other.dimension_ && moduli_ == other.moduli_ &&
               components_ == other.components_;
    }

    /**
     * @brief Returns the number of components, i.e. the number of moduli in the
     * coefficient RNS.
     * @return const size_t
     */
    inline const size_t component_count() const { return components_.size(); }

    /**
     * @brief Returns log2 value of component polynomials' length.
     * @return const size_t
     */
    inline const size_t log_dimension() const { return log_dimension_; }

    /**
     * @brief Returns component polynomials' length.
     * @return const size_t
     */
    inline const size_t dimension() const { return dimension_; }

    /**
     * @brief Returns the components_ vector as non-const reference.
     * @return std::vector<ComponentData> &
     */
    inline std::vector<ComponentData> &components() { return components_; }

    /**
     * @brief Returns the components_ vector as const reference.
     * @return const std::vector<ComponentData> &
     */
    inline const std::vector<ComponentData> &components() const {
        return components_;
    }

    /**
     * @brief Returns the begin of components_ vector as non-const iterator.
     * @return std::vector<ComponentData>::iterator
     */
    inline auto begin() { return components_.begin(); }

    /**
     * @brief Returns the begin of components_ vector as const iterator.
     * @return const std::vector<ComponentData>::const_iterator
     */
    inline const auto begin() const { return components_.cbegin(); }

    /**
     * @brief Returns the end of components_ vector as non-const iterator.
     * @return std::vector<ComponentData>::iterator
     */
    inline auto end() { return components_.end(); }

    /**
     * @brief Returns the end of components_ vector as const iterator.
     * @return const std::vector<ComponentData>::const_iterator
     */
    inline const auto end() const { return components_.cend(); }

    /**
     * @brief Returns the last item of components_ vector as non-const iterator.
     * @return std::vector<ComponentData>::iterator
     */
    inline auto last() { return components_.end() - 1; }

    /**
     * @brief Returns the last item of components_ vector as const iterator.
     * @return const std::vector<ComponentData>::const_iterator
     */
    inline const auto last() const { return components_.cend() - 1; }

    /**
     * @brief Get the i-th modulus in moduli_.
     * @param i The index of the accessing modulus.
     * @return const u64
     */
    inline const u64 modulus_at(int i) const { return moduli_[i]; }

    /**
     * @brief Get the moduli_ vector.
     * @return const std::vector<u64> &
     */
    inline const std::vector<u64> &modulus_vec() const { return moduli_; }

    /**
     * @brief Get the i-th component as non-const reference.
     * @param i The index of the accessing component.
     * @return u64 *
     */
    inline ComponentData &operator[](int i) { return components_[i]; }

    /**
     * @brief Get the i-th component as const reference.
     * @param i The index of the accessing component.
     * @return const u64 *
     */
    inline const ComponentData &operator[](int i) const {
        return components_[i];
    }

    /**
     * @brief Modify the RnsIntVec to enlarge the RNS, i.e. add some
     * components and corresponding moduli.
     * @param new_moduli The new moduli to be added into the RNS basis.
     * @param adding The number of components which are to be added.
     */
    void add_components(const std::vector<u64> &new_moduli, size_t adding = 1);

    /**
     * @brief Modify the RnsIntVec to shrink the RNS, i.e. remove some
     * components and corresponding moduli.
     * @param removing The number of components which are to be removed.
     */
    void remove_components(size_t removing = 1);

    friend const RnsIntVec &operator+=(RnsIntVec &self, const RnsIntVec &b);

    friend RnsIntVec operator+(const RnsIntVec &a, const RnsIntVec &b);

    friend const RnsIntVec &operator-=(RnsIntVec &self, const RnsIntVec &b);

    friend RnsIntVec operator-(const RnsIntVec &a, const RnsIntVec &b);

    friend const RnsIntVec &operator*=(RnsIntVec &self, const RnsIntVec &b);

    friend RnsIntVec operator*(const RnsIntVec &a, const RnsIntVec &b);

    friend const RnsIntVec &operator*=(RnsIntVec &self, const u64 small_scalar);

    friend const RnsIntVec &operator*=(RnsIntVec &self,
                                       const std::vector<u64> &rns_scalar);

private:
    /// log2 value of component polynomials' length.
    size_t log_dimension_ = 0;

    /// Component polynomials' length.
    size_t dimension_ = 0;

    /// A vector storing all polynomial coefficients.
    std::vector<ComponentData> components_;

    /// A vector of RNS basis, i.e. all the moduli of coefficients.
    std::vector<u64> moduli_;
};

class RnsPolynomial : public RnsIntVec {
public:
    using RnsIntVec::RnsIntVec;

    /**
     * @brief Representation form of the polynomial.
     */
    enum class RepForm { coeff, value };

    RnsPolynomial(RnsIntVec &&rns_int_vec) : RnsIntVec(rns_int_vec) {}

    friend void ntt_negacyclic_inplace_lazy(RnsPolynomial &);

    friend void intt_negacyclic_inplace_lazy(RnsPolynomial &);

    friend const RnsPolynomial &operator+=(RnsPolynomial &self,
                                           const RnsPolynomial &b);

    friend RnsPolynomial operator+(const RnsPolynomial &a, const RnsPolynomial &b);

    friend const RnsPolynomial &operator-=(RnsPolynomial &self,
                                           const RnsPolynomial &b);

    friend RnsPolynomial operator-(const RnsPolynomial &a, const RnsPolynomial &b);

    friend const RnsPolynomial &operator*=(RnsPolynomial &self,
                                           const RnsPolynomial &b);

    friend RnsPolynomial operator*(const RnsPolynomial &a,
                                   const RnsPolynomial &b);

    friend const RnsPolynomial &operator*=(RnsPolynomial &self,
                                           const u64 small_scalar);

    friend const RnsPolynomial &operator*=(RnsPolynomial &self,
                                           const std::vector<u64> &rns_scalar);

    /// Representation form of the polynomial, default being coefficients.
    /// This is set to be publicly visible in order to enable possible tweaks.
    RepForm rep_form = RepForm::coeff;
};

using RnsPolyParams = RnsPolynomial::Params;

using PolyRepForm = RnsPolynomial::RepForm;

/**
 * @brief TODO
 *
 * @param self
 * @param b
 * @return RnsIntVec
 */
const RnsIntVec &operator+=(RnsIntVec &self, const RnsIntVec &b);

/**
 * @brief TODO
 *
 * @param a
 * @param b
 * @return RnsIntVec
 */
inline RnsIntVec operator+(const RnsIntVec &a, const RnsIntVec &b) {
    auto result(a);
    result += b;
    return result;
}

/**
 * @brief TODO
 *
 * @param self
 * @param b
 * @return RnsIntVec
 */
const RnsIntVec &operator-=(RnsIntVec &self, const RnsIntVec &b);

/**
 * @brief TODO
 *
 * @param a
 * @param b
 * @return RnsIntVec
 */
inline RnsIntVec operator-(const RnsIntVec &a, const RnsIntVec &b) {
    auto result(a);
    result -= b;
    return result;
}

/**
 * @brief TODO
 *
 * @param a
 * @param b
 * @return RnsIntVec
 */
RnsIntVec operator*(const RnsIntVec &a, const RnsIntVec &b);

/**
 * @brief TODO
 *
 * @param self
 * @param b
 * @return RnsIntVec
 */
inline const RnsIntVec &operator*=(RnsIntVec &self, const RnsIntVec &b) {
    auto temp(self);
    return self = temp * b;
}

/**
 * @brief TODO
 *
 * @param self
 * @param small_scalar
 * @return const RnsIntVec&
 */
const RnsIntVec &operator*=(RnsIntVec &self, const u64 small_scalar);

/**
 * @brief TODO
 *
 * @param poly
 * @param small_scalar
 * @return RnsIntVec&
 */
inline RnsIntVec operator*(const RnsIntVec &int_vec, const u64 small_scalar) {
    auto int_vec_copy(int_vec);
    int_vec_copy *= small_scalar;
    return int_vec_copy;
}

/**
 * @brief TODO
 *
 * @param self
 * @param rns_scalar
 * @return const RnsIntVec&
 */
const RnsIntVec &operator*=(RnsIntVec &self,
                            const std::vector<u64> &rns_scalar);

/**
 * @brief TODO
 *
 * @param poly
 * @param rns_scalar
 * @return RnsIntVec
 */
inline RnsIntVec operator*(const RnsIntVec &int_vec,
                           const std::vector<u64> &rns_scalar) {
    auto int_vec_copy(int_vec);
    int_vec_copy *= rns_scalar;
    return int_vec_copy;
}

#ifdef HEHUB_DEBUG_FHE
/**
 * @brief TODO
 *
 * @param out
 * @param rns_poly
 * @return std::ostream&
 */
std::ostream &operator<<(std::ostream &out, const RnsIntVec &rns_poly);
#endif

inline const RnsPolynomial &operator+=(RnsPolynomial &self,
                                       const RnsPolynomial &b) {
    if (self.rep_form != b.rep_form) {
        throw std::invalid_argument(
            "Operands are in different representation form.");
    }

    self += (RnsIntVec &)b;
    return self;
}

inline RnsPolynomial operator+(const RnsPolynomial &a, const RnsPolynomial &b) {
    auto result(a);
    result += b;
    return result;
}

inline const RnsPolynomial &operator-=(RnsPolynomial &self,
                                       const RnsPolynomial &b) {
    if (self.rep_form != b.rep_form) {
        throw std::invalid_argument(
            "Operands are in different representation form.");
    }

    self -= (const RnsIntVec &)b;
    return self;
}

inline RnsPolynomial operator-(const RnsPolynomial &a, const RnsPolynomial &b) {
    auto result(a);
    result -= b;
    return result;
}

inline RnsPolynomial operator*(const RnsPolynomial &a, const RnsPolynomial &b) {
    if (a.rep_form == PolyRepForm::coeff) {
        throw std::invalid_argument("Operand a is in coefficient form.");
    }
    if (b.rep_form == PolyRepForm::coeff) {
        throw std::invalid_argument("Operand b is in coefficient form.");
    }

    RnsPolynomial result = (const RnsIntVec &)a * (const RnsIntVec &)b;
    result.rep_form = PolyRepForm::value;

    return result;
}

inline const RnsPolynomial &operator*=(RnsPolynomial &self,
                                       const RnsPolynomial &b) {
    auto temp(self);
    return self = temp * b;
}

inline const RnsPolynomial &operator*=(RnsPolynomial &self,
                                       const u64 small_scalar) {
    RnsIntVec &self_ref = self;
    self_ref *= small_scalar;
    return self;
}

inline const RnsPolynomial &operator*=(RnsPolynomial &self,
                                       const std::vector<u64> &rns_scalar) {
    RnsIntVec &self_ref = self;
    self_ref *= rns_scalar;
    return self;
}

inline RnsPolynomial operator*(const RnsPolynomial &poly,
                               const std::vector<u64> &rns_scalar) {
    auto poly_copy(poly);
    poly_copy *= rns_scalar;
    return poly_copy;
}

} // namespace hehub
