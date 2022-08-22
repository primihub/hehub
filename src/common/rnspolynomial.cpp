#include "rnspolynomial.h"
#include <iostream>

namespace hehub {

RnsPolynomial::ComponentData::ComponentData(const size_t poly_len)
    : poly_len_(poly_len), data_(new u64[poly_len]) {}

RnsPolynomial::ComponentData::ComponentData(const ComponentData &other)
    : poly_len_(other.poly_len_), data_(new u64[other.poly_len_]) {
    std::copy(other.data_, other.data_ + poly_len_, data_);
}

RnsPolynomial::ComponentData::~ComponentData() {
    if (data_ != nullptr) {
        delete[] data_;
    }
}

RnsPolynomial::ComponentData &
RnsPolynomial::ComponentData::operator=(const ComponentData &copying) {
    if (poly_len_ != copying.poly_len_) {
        if (data_ != nullptr) {
            delete[] data_;
        }
        poly_len_ = copying.poly_len_;
        data_ = new u64[copying.poly_len_];
    }
    std::copy(copying.data_, copying.data_ + poly_len_, data_);
    return *this;
}

RnsPolynomial::ComponentData &
RnsPolynomial::ComponentData::operator=(ComponentData &&moving) noexcept {
    if (this == &moving) {
        return *this;
    }

    poly_len_ = moving.poly_len_;
    if (data_ != nullptr) {
        delete[] data_;
    }
    data_ = moving.data_;
    moving.data_ = nullptr;
    return *this;
}

RnsPolynomial::RnsPolynomial(const size_t components, const size_t log_poly_len,
                             const std::vector<u64> &moduli)
    : components_(components), log_poly_len_(log_poly_len),
      poly_len_(1 << log_poly_len) {

    if (moduli.size() < component_count()) {
        throw std::invalid_argument(
            "No matching number of moduli provided to create RnsPolynomial.");
    }
    moduli_.assign(moduli.begin(), moduli.begin() + component_count());
    for (auto &component: components_) {
        component = ComponentData(poly_len_);
    }
}

void RnsPolynomial::add_components(const std::vector<u64> &new_moduli,
                                   size_t adding) {
    if (new_moduli.size() < adding) {
        throw std::invalid_argument(
            "No matching number of moduli provided to add components.");
    }

    moduli_.insert(moduli_.end(), new_moduli.begin(), new_moduli.end());
    components_.resize(components_.size() + adding);
}

void RnsPolynomial::remove_components(size_t removing) {
    if (component_count() < removing) {
        throw std::invalid_argument(
            "Trying to remove components more than existing.");
    }

    moduli_.erase(moduli_.end() - removing, moduli_.end());
    components_.erase(components_.end() - removing, components_.end());
}

} // namespace hehub
