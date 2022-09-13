#include "rnspolynomial.h"
#include "mod_arith.h"
#include <cmath>

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

RnsPolynomial::RnsPolynomial(const size_t poly_len, const size_t components,
                             const std::vector<u64> &moduli)
    : poly_len_(poly_len), components_(components),
      log_poly_len_(std::log2(poly_len) + 0.5) {

    if (poly_len_ != 1 << log_poly_len_) {
        throw std::invalid_argument("poly_len should be a 2-power.");
    }

    if (moduli.size() < component_count()) {
        throw std::invalid_argument(
            "No matching number of moduli provided to create RnsPolynomial.");
    }
    moduli_.assign(moduli.begin(), moduli.begin() + component_count());
    for (auto &component : components_) {
        component = ComponentData(poly_len_);
    }
}

RnsPolynomial::RnsPolynomial(const PolyDimensions &poly_dim)
    : RnsPolynomial(poly_dim.poly_len, poly_dim.component_count,
                    poly_dim.moduli) {}

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

const RnsPolynomial &operator+=(RnsPolynomial &self, const RnsPolynomial &b) {
    if (self.rep_form != b.rep_form) {
        throw std::invalid_argument(
            "Operands are in different representation form.");
    }
    if (self.poly_len() != b.poly_len()) {
        throw std::invalid_argument("Operands' poly len mismatch.");
    }
    auto poly_len = self.poly_len();
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
        for (size_t i = 0; i < poly_len; i++) {
            self[k][i] += b[k][i];
            self[k][i] -=
                (self[k][i] >= moduli_doubled[k]) ? moduli_doubled[k] : 0;
        }
    }

    return self;
}

const RnsPolynomial &operator-=(RnsPolynomial &self, const RnsPolynomial &b) {
    if (self.rep_form != b.rep_form) {
        throw std::invalid_argument(
            "Operands are in different representation form.");
    }
    if (self.poly_len() != b.poly_len()) {
        throw std::invalid_argument("Operands' poly len mismatch.");
    }
    auto poly_len = self.poly_len();
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
        for (size_t i = 0; i < poly_len; i++) {
            self[k][i] += moduli_doubled[k] - b[k][i];
            self[k][i] -=
                (self[k][i] >= moduli_doubled[k]) ? moduli_doubled[k] : 0;
        }
    }

    return self;
}

RnsPolynomial operator*(const RnsPolynomial &a, const RnsPolynomial &b) {
    if (a.rep_form == PolyRepForm::coeff) {
        throw std::invalid_argument("Operand a is in coefficient form.");
    }
    if (b.rep_form == PolyRepForm::coeff) {
        throw std::invalid_argument("Operand b is in coefficient form.");
    }
    if (a.poly_len() != b.poly_len()) {
        throw std::invalid_argument("Operands' poly len mismatch.");
    }
    auto poly_len = a.poly_len();
    auto components = std::min(a.component_count(), b.component_count());
    auto moduli(a.modulus_vec()), b_moduli(b.modulus_vec());
    moduli.resize(components);
    b_moduli.resize(components);
    if (moduli != b_moduli) {
        throw std::invalid_argument("Operands' moduli mismatch.");
    }

    RnsPolynomial result(PolyDimensions{poly_len, components, moduli});
    result.rep_form = PolyRepForm::value;
    for (size_t k = 0; k < components; k++) {
        batched_mul_mod_hybrid_lazy(moduli[k], poly_len, a[k].data(), b[k].data(), result[k].data());
    }

    return result;
}

#ifdef FHE_DEBUG
std::ostream &operator<<(std::ostream &out, const RnsPolynomial &rns_poly) {
    auto component_count = rns_poly.component_count();
    auto poly_len = rns_poly.poly_len();
    auto mod_ptr = rns_poly.modulus_vec().begin();

    for (const auto &component_poly : rns_poly) {
        out << "mod " << *(mod_ptr++) << ":\t[ ";
        for (const auto &coeff : component_poly) {
            // here "coeff" can also mean NTT value
            out << coeff << ", ";
        }
        out << "]" << std::endl;
    }

    return out;
}
#endif

} // namespace hehub
