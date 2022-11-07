/**
 * @file rns.h
 * @brief TODO
 *
 */
#pragma once

#include "cereal/archives/binary.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/vector.hpp"
#include "type_defs.h"
#include <iomanip>
#include <sstream>
#include <vector>

namespace hehub {

class RnsIntVec {
public:
    struct Params {
        size_t dimension = 0;

        size_t component_count;

        std::vector<u64> moduli;
    };

    class ComponentData {
    public:
        ComponentData() {}

        ComponentData(const size_t dimension);

        ComponentData(const ComponentData &other);

        ComponentData(ComponentData &&other) noexcept {
            *this = std::move(other);
        }

        ~ComponentData();

        inline u64 &operator[](const int idx) { return data_[idx]; }

        inline const u64 operator[](const int idx) const { return data_[idx]; }

        ComponentData &operator=(const ComponentData &copying);

        ComponentData &operator=(ComponentData &&moving) noexcept;

        bool operator==(const ComponentData &other) const;

        inline bool operator!=(const ComponentData &comparing) const {
            return !((*this) == comparing);
        }

        inline u64 *data() { return data_; }

        inline const u64 *data() const { return data_; }

        inline u64 *begin() { return data_; }

        inline const u64 *begin() const { return data_; }

        inline u64 *end() { return data_ + dimension_; }

        inline const u64 *end() const { return data_ + dimension_; }

        template <class Archive> std::string save_minimal(Archive &) const {
            std::stringstream ss;
            ss << std::hex << std::setw(8) << dimension_;
            std::string archived = ss.str();
            archived.resize(dimension_ * 8 + 8, 0);
            std::copy((char *)data_, (char *)data_ + dimension_ * 8, archived.begin() + 8);
            return archived;
        }

        template <class Archive>
        void load_minimal(const Archive &, const std::string &archived) {
            if (data_ != nullptr) {
                delete [] data_;
            }
            dimension_ = std::stoi(archived.substr(0, 8), nullptr, 16);
            if (archived.size() != dimension_ * 8 + 8) {
                throw std::invalid_argument("Archived data format invalid.");
            }
            data_ = new u64[dimension_];
            std::copy(archived.begin() + 8, archived.end(), (char *)data_);
        }

    private:
        u64 *data_ = nullptr;

        size_t dimension_ = 0;
    };

    enum class RepForm { coeff, value };

    RnsIntVec() {}

    RnsIntVec(const size_t dimension, const size_t components,
              const std::vector<u64> &moduli);

    RnsIntVec(const Params &params);

    RnsIntVec(const RnsIntVec &other) = default;

    RnsIntVec(RnsIntVec &&other) = default;

    RnsIntVec &operator=(const RnsIntVec &other) = default;

    RnsIntVec &operator=(RnsIntVec &&other) = default;

    inline const bool operator==(const RnsIntVec &other) const {
        return log_dimension_ == other.log_dimension_ &&
               dimension_ == other.dimension_ && moduli_ == other.moduli_ &&
               components_ == other.components_;
    }

    inline const size_t component_count() const { return components_.size(); }

    inline const size_t log_dimension() const { return log_dimension_; }

    inline const size_t dimension() const { return dimension_; }

    inline std::vector<ComponentData> &components() { return components_; }

    inline const std::vector<ComponentData> &components() const {
        return components_;
    }

    inline auto begin() { return components_.begin(); }

    inline const auto begin() const { return components_.cbegin(); }

    inline auto end() { return components_.end(); }

    inline const auto end() const { return components_.cend(); }

    inline auto last() { return components_.end() - 1; }

    inline const auto last() const { return components_.cend() - 1; }

    inline const u64 modulus_at(int i) const { return moduli_[i]; }

    inline const std::vector<u64> &modulus_vec() const { return moduli_; }

    inline ComponentData &operator[](int i) { return components_[i]; }

    inline const ComponentData &operator[](int i) const {
        return components_[i];
    }

    void add_components(const std::vector<u64> &new_moduli, size_t adding = 1);

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

    template <class Archive> void save(Archive &ar) const {
        ar(log_dimension_, dimension_, components_, moduli_);
    }

    template <class Archive> void load(Archive &ar) {
        ar(log_dimension_, dimension_, components_, moduli_);
    }

private:
    size_t log_dimension_ = 0;

    size_t dimension_ = 0;

    std::vector<ComponentData> components_;

    std::vector<u64> moduli_;
};

class RnsPolynomial : public RnsIntVec {
public:
    using RnsIntVec::RnsIntVec;

    enum class RepForm { coeff, value };

    RnsPolynomial(RnsIntVec &&rns_int_vec) : RnsIntVec(rns_int_vec) {}

    friend void ntt_negacyclic_inplace_lazy(RnsPolynomial &);

    friend void intt_negacyclic_inplace_lazy(RnsPolynomial &);

    friend const RnsPolynomial &operator+=(RnsPolynomial &self,
                                           const RnsPolynomial &b);

    friend RnsPolynomial operator+(const RnsPolynomial &a,
                                   const RnsPolynomial &b);

    friend const RnsPolynomial &operator-=(RnsPolynomial &self,
                                           const RnsPolynomial &b);

    friend RnsPolynomial operator-(const RnsPolynomial &a,
                                   const RnsPolynomial &b);

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

const RnsIntVec &operator+=(RnsIntVec &self, const RnsIntVec &b);

inline RnsIntVec operator+(const RnsIntVec &a, const RnsIntVec &b) {
    auto result(a);
    result += b;
    return result;
}

const RnsIntVec &operator-=(RnsIntVec &self, const RnsIntVec &b);

inline RnsIntVec operator-(const RnsIntVec &a, const RnsIntVec &b) {
    auto result(a);
    result -= b;
    return result;
}

RnsIntVec operator*(const RnsIntVec &a, const RnsIntVec &b);

inline const RnsIntVec &operator*=(RnsIntVec &self, const RnsIntVec &b) {
    auto temp(self);
    return self = temp * b;
}

const RnsIntVec &operator*=(RnsIntVec &self, const u64 small_scalar);

inline RnsIntVec operator*(const RnsIntVec &int_vec, const u64 small_scalar) {
    auto int_vec_copy(int_vec);
    int_vec_copy *= small_scalar;
    return int_vec_copy;
}

const RnsIntVec &operator*=(RnsIntVec &self,
                            const std::vector<u64> &rns_scalar);

inline RnsIntVec operator*(const RnsIntVec &int_vec,
                           const std::vector<u64> &rns_scalar) {
    auto int_vec_copy(int_vec);
    int_vec_copy *= rns_scalar;
    return int_vec_copy;
}

#ifdef HEHUB_DEBUG_FHE
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
