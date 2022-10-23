#include "rns.h"
#include "mod_arith.h"
#include "range/v3/view/zip.hpp"
#include <cmath>

using namespace ranges::views;

namespace hehub {

RnsIntVec::ComponentData::ComponentData(const size_t dimension)
    : dimension_(dimension), data_(new u64[dimension]) {}

RnsIntVec::ComponentData::ComponentData(const ComponentData &other)
    : dimension_(other.dimension_), data_(new u64[other.dimension_]) {
    std::copy(other.data_, other.data_ + dimension_, data_);
}

RnsIntVec::ComponentData::~ComponentData() {
    if (data_ != nullptr) {
        delete[] data_;
    }
}

RnsIntVec::ComponentData &
RnsIntVec::ComponentData::operator=(const ComponentData &copying) {
    if (dimension_ != copying.dimension_) {
        if (data_ != nullptr) {
            delete[] data_;
        }
        dimension_ = copying.dimension_;
        data_ = new u64[copying.dimension_];
    }
    std::copy(copying.data_, copying.data_ + dimension_, data_);
    return *this;
}

RnsIntVec::ComponentData &
RnsIntVec::ComponentData::operator=(ComponentData &&moving) noexcept {
    if (this == &moving) {
        return *this;
    }

    dimension_ = moving.dimension_;
    if (data_ != nullptr) {
        delete[] data_;
    }
    data_ = moving.data_;
    moving.data_ = nullptr;
    return *this;
}

RnsIntVec::RnsIntVec(const size_t dimension, const size_t components,
                     const std::vector<u64> &moduli)
    : dimension_(dimension), components_(components),
      log_dimension_(std::log2(dimension) + 0.5) {

    // This condition should be moved to RnsPolynomial
    if (dimension_ != 1 << log_dimension_) {
        throw std::invalid_argument("dimension should be a 2-power.");
    }

    if (moduli.size() < component_count()) {
        throw std::invalid_argument(
            "No matching number of moduli provided to create RnsIntVec.");
    }
    moduli_.assign(moduli.begin(), moduli.begin() + component_count());
    for (auto &component : components_) {
        component = ComponentData(dimension_);
    }
}

RnsIntVec::RnsIntVec(const RnsIntVec::Params &params)
    : RnsIntVec(params.dimension, params.component_count, params.moduli) {}

void RnsIntVec::add_components(const std::vector<u64> &new_moduli,
                               size_t adding) {
    if (new_moduli.size() < adding) {
        throw std::invalid_argument(
            "No matching number of moduli provided to add components.");
    }

    auto orig_size = components_.size();
    moduli_.insert(moduli_.end(), new_moduli.begin(), new_moduli.end());
    components_.resize(components_.size() + adding);
    for (size_t i = orig_size; i < components_.size(); i++) {
        components_[i] = ComponentData(dimension_);
    }
}

void RnsIntVec::remove_components(size_t removing) {
    if (component_count() < removing) {
        throw std::invalid_argument(
            "Trying to remove components more than existing.");
    }

    moduli_.erase(moduli_.end() - removing, moduli_.end());
    components_.erase(components_.end() - removing, components_.end());
}

const RnsIntVec &operator+=(RnsIntVec &self, const RnsIntVec &b) {
    if (self.dimension() != b.dimension()) {
        throw std::invalid_argument("Operands' poly len mismatch.");
    }
    auto dimension = self.dimension();
    if (b.component_count() < self.component_count()) {
        throw std::invalid_argument(
            "Operand b contains less components than self.");
    }
    auto components = self.component_count();
    auto moduli(self.modulus_vec()), b_moduli(b.modulus_vec());
    b_moduli.resize(components);
    if (moduli != b_moduli) {
        throw std::invalid_argument("Operands' moduli mismatch.");
    }

    auto moduli_doubled(std::move(moduli));
    for (auto &m : moduli_doubled) {
        m *= 2;
    }
    for (size_t k = 0; k < components; k++) {
        for (size_t i = 0; i < dimension; i++) {
            self[k][i] += b[k][i];
            self[k][i] -=
                (self[k][i] >= moduli_doubled[k]) ? moduli_doubled[k] : 0;
        }
    }

    return self;
}

const RnsIntVec &operator-=(RnsIntVec &self, const RnsIntVec &b) {
    if (self.dimension() != b.dimension()) {
        throw std::invalid_argument("Operands' poly len mismatch.");
    }
    auto dimension = self.dimension();
    if (b.component_count() < self.component_count()) {
        throw std::invalid_argument(
            "Operand b contains less components than self.");
    }
    auto components = self.component_count();
    auto moduli(self.modulus_vec()), b_moduli(b.modulus_vec());
    b_moduli.resize(components);
    if (moduli != b_moduli) {
        throw std::invalid_argument("Operands' moduli mismatch.");
    }

    auto moduli_doubled(std::move(moduli));
    for (auto &m : moduli_doubled) {
        m *= 2;
    }
    for (size_t k = 0; k < components; k++) {
        for (size_t i = 0; i < dimension; i++) {
            self[k][i] += moduli_doubled[k] - b[k][i];
            self[k][i] -=
                (self[k][i] >= moduli_doubled[k]) ? moduli_doubled[k] : 0;
        }
    }

    return self;
}

RnsIntVec operator*(const RnsIntVec &a, const RnsIntVec &b) {
    if (a.dimension() != b.dimension()) {
        throw std::invalid_argument("Operands' poly len mismatch.");
    }
    auto dimension = a.dimension();
    auto components = std::min(a.component_count(), b.component_count());
    auto moduli(a.modulus_vec()), b_moduli(b.modulus_vec());
    moduli.resize(components);
    b_moduli.resize(components);
    if (moduli != b_moduli) {
        throw std::invalid_argument("Operands' moduli mismatch.");
    }

    RnsIntVec result(RnsIntVec::Params{dimension, components, moduli});
    for (size_t k = 0; k < components; k++) {
        batched_mul_mod_hybrid_lazy(moduli[k], dimension, a[k].data(),
                                    b[k].data(), result[k].data());
    }

    return result;
}

const RnsIntVec &operator*=(RnsIntVec &self, const u64 small_scalar) {
    for (size_t k = 0; k < self.component_count(); k++) {
        auto curr_mod = self.moduli_[k];
        auto scalar_reduced = small_scalar % curr_mod; // need opt?
        auto scalar_harvey = ((u128)scalar_reduced << 64) / curr_mod;
        for (auto &coeff : self[k]) {
            coeff = mul_mod_harvey_lazy(curr_mod, coeff, scalar_reduced,
                                        scalar_harvey);
        }
    }
    return self;
}

const RnsIntVec &operator*=(RnsIntVec &self,
                            const std::vector<u64> &rns_scalar) {
    if (rns_scalar.size() != self.component_count()) {
        throw std::invalid_argument("Numbers of RNS component mismatch.");
    }

    for (size_t k = 0; k < self.component_count(); k++) {
        auto curr_mod = self.moduli_[k];
        auto scalar_reduced = rns_scalar[k] % curr_mod; // need opt?
        auto scalar_harvey = ((u128)scalar_reduced << 64) / curr_mod;
        for (auto &coeff : self[k]) {
            coeff = mul_mod_harvey_lazy(curr_mod, coeff, scalar_reduced,
                                        scalar_harvey);
        }
    }
    return self;
}

#ifdef HEHUB_DEBUG_FHE
std::ostream &operator<<(std::ostream &out, const RnsIntVec &rns_poly) {
    auto component_count = rns_poly.component_count();
    auto dimension = rns_poly.dimension();
    auto &moduli = rns_poly.modulus_vec();

    for (auto [component, modulus] : zip(rns_poly, moduli)) {
        out << "mod " << modulus << ":\t[ ";
        for (const auto &coeff : component) {
            // here "coeff" can also mean NTT value
            out << coeff << ", ";
        }
        out << "]" << std::endl;
    }

    return out;
}
#endif

} // namespace hehub
