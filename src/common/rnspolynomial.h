/**
 * @file RnsPolynomial.h
 * @brief RnsPolynomial class.
 *
 */
#pragma once

#include "type_defs.h"
#include <sstream>
#include <vector>

namespace hehub {

/**
 * @class RnsPolynomial
 * @brief A RnsPolynomial mainly contains a group of polynomials whose
 * coefficients of terms of same degree constitute an residue number system
 * (RNS). The degree of polynomials equals certain 2-power minus 1.
 */
class RnsPolynomial {
public:
    /**
     * @class ComponentData
     * @brief A dynamic array for storing all coefficients of a component
     * polynomial of a RnsPolynomial object.
     */
    class ComponentData {
    public:
        /**
         * @brief Construct a new Simple Poly object.
         */
        ComponentData() {}

        /**
         * @brief Construct a new ComponentData object.
         * @param poly_len The length of the polynomial.
         */
        ComponentData(const size_t poly_len);

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
            if (poly_len_ != other.poly_len_)
                return false;
            for (int i = 0; i < poly_len_; i++) {
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

    private:
        /// Pointer as a dynamic array for storing the coefficients.
        u64 *data_ = nullptr;

        /// The length of the polynomial (= degree + 1).
        size_t poly_len_ = 0;
    };

    /**
     * @brief Construct an empty RnsPolynomial object.
     */
    RnsPolynomial() {}

    /**
     * @brief Construct a new RnsPolynomial object.
     * @param components The number of components, i.e. polynomials each
     * modulo different integer.
     * @param log_poly_len log2 value of each component polynomial's length
     * (= degree + 1).
     */
    RnsPolynomial(const size_t components, const size_t log_poly_len,
                  const std::vector<u64> &moduli);

    /**
     * @brief Creates an RnsPolynomial object by copying from another one.
     * @param other The source RnsPolynomial object.
     */
    RnsPolynomial(const RnsPolynomial &other) = default;

    /**
     * @brief Creates an RnsPolynomial object by moving another one.
     * @param other The source RnsPolynomial object.
     */
    RnsPolynomial(RnsPolynomial &&other) = default;

    /**
     * @brief Copies a given RnsPolynomial to this.
     * @param other The given RnsPolynomial object.
     * @return RnsPolynomial &
     */
    RnsPolynomial &operator=(const RnsPolynomial &other) = default;

    /**
     * @brief Moves a given RnsPolynomial to replace this.
     * @param other The given RnsPolynomial object.
     * @return RnsPolynomial &
     */
    RnsPolynomial &operator=(RnsPolynomial &&other) = default;

    /**
     * @brief Compares this with another RnsPolynomial, returning true only when
     * all spec and data coincide.
     * @param other The RnsPolynomial object for comparison.
     * @return const bool
     */
    inline const bool operator==(const RnsPolynomial &other) const {
        return log_poly_len_ == other.log_poly_len_ &&
               poly_len_ == other.poly_len_ && moduli_ == other.moduli_ &&
               components_ == other.components_;
    }

    /**
     * @brief Returns the number of components, i.e. the number of moduli in the
     * coefficient RNS.
     * @return const size_t
     */
    inline const size_t component_count() const { return components_.size(); }

    /**
     * @brief Returns log2 value of component polynomials' length (which is
     * their degree + 1).
     * @return const size_t
     */
    inline const size_t log_poly_len() const { return log_poly_len_; }

    /**
     * @brief Returns component polynomials' length (= their
     * degree + 1).
     * @return const size_t
     */
    inline const size_t poly_len() const { return poly_len_; }

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
     * @brief Get the i-th modulus in moduli_.
     * @param i The index of the accessing modulus.
     * @return const u64
     */
    inline const u64 modulus_at(int i) const { return moduli_[i]; }

    /**
     * @brief Get the moduli_ vector.
     * @return const std::vector<u64> &
     */
    inline const std::vector<u64> &moduli_vec() const { return moduli_; }

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
     * @brief Modify the RnsPolynomial to enlarge the RNS, i.e. add some
     * components and corresponding moduli.
     * @param new_moduli The new moduli to be added into the RNS basis.
     * @param adding The number of components which are to be added.
     */
    void add_components(const std::vector<u64> &new_moduli, size_t adding = 1);

    /**
     * @brief Modify the RnsPolynomial to shrink the RNS, i.e. remove some
     * components and corresponding moduli.
     * @param removing The number of components which are to be removed.
     */
    void remove_components(size_t removing = 1);

    // void save(std::stringstream &stream);

    // void load(std::stringstream &stream, u64 log_poly_len);

private:
    /// log2 value of component polynomials' length (which is their degree + 1).
    size_t log_poly_len_ = 0;

    /// Component polynomials' length (= their degree + 1).
    size_t poly_len_ = 0;

    /// A vector storing all polynomial coefficients.
    std::vector<ComponentData> components_;

    /// A vector of RNS basis, i.e. all the moduli of coefficients.
    std::vector<u64> moduli_;
};

} // namespace hehub
