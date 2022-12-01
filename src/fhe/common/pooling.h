#pragma once

#include "type_defs.h"
#include <cassert>
#include <map>
#include <set>
#include <stack>
#include <type_traits>

namespace hehub {

class Poolable {
public:
    Poolable() {}

    virtual void release() = 0;
};

// template <typename T, typename std::enable_if_t<
//                           std::is_pointer_v<T> &&
//                           std::is_base_of_v<Poolable,
//                           std::remove_pointer_t<T>>>
//                           * = nullptr>
template <typename T,
          typename std::enable_if_t<std::is_pointer_v<T>> * = nullptr>
struct SimpleObjPool {
    std::stack<T> cached_list;

    std::set<T> objects;

    ~SimpleObjPool() {
        std::stack<T>().swap(cached_list);

        for (auto obj_ptr : objects) {
            obj_ptr->release();
        }
    }
};

template <typename T> class AutopoolArray : public Poolable {
    static std::map<size_t, SimpleObjPool<AutopoolArray<T> *>> general_pool_;

public:
    AutopoolArray() {}

    AutopoolArray(const size_t dimension) { require(dimension); }

    AutopoolArray(const AutopoolArray &other) {
        require(other.dimension_);
        std::copy(other.data_, other.data_ + dimension_, data_);
    }

    AutopoolArray(AutopoolArray &&other) noexcept { *this = std::move(other); }

    ~AutopoolArray() { cache(); }

    inline T &operator[](const int idx) { return data_[idx]; }

    inline const T operator[](const int idx) const { return data_[idx]; }

    AutopoolArray &operator=(const AutopoolArray &copying) {
        if (dimension_ != copying.dimension_) {
            cache();
            require(copying.dimension_);
        }
        std::copy(copying.data_, copying.data_ + dimension_, data_);
        return *this;
    }

    AutopoolArray &operator=(AutopoolArray &&moving) noexcept {
        if (this == &moving) {
            return *this;
        }

        aff_pool_ = moving.aff_pool_;
        if (!aff_pool_) {
            init_aff_pool(moving.dimension_);
        }
        if (aff_pool_->objects.find(this) != aff_pool_->objects.end()) {
            cache();
        } else {
            aff_pool_->objects.insert(this);
        }

        dimension_ = moving.dimension_;
        moving.dimension_ = 0;
        data_ = moving.data_;
        moving.data_ = nullptr;
        assert(aff_pool_->objects.find(&moving) != aff_pool_->objects.end());
        aff_pool_->objects.erase(aff_pool_->objects.find(&moving));

        return *this;
    }

    inline bool operator==(const AutopoolArray &other) const {
        if (dimension_ != other.dimension_)
            return false;
        for (int i = 0; i < dimension_; i++) {
            if ((*this)[i] != other[i]) {
                return false;
            }
        }
        return true;
    }

    inline bool operator!=(const AutopoolArray &comparing) const {
        return !((*this) == comparing);
    }

    inline T *data() { return data_; }

    inline const T *data() const { return data_; }

    inline T *begin() { return data_; }

    inline const T *begin() const { return data_; }

    inline T *end() { return data_ + dimension_; }

    inline const T *end() const { return data_ + dimension_; }

    inline const auto &aff_pool() const { return *aff_pool_; }

    void init_aff_pool(size_t dimension) {
        auto pool_iter = general_pool_.find(dimension);
        if (pool_iter == general_pool_.end()) {
            general_pool_.insert(std::make_pair(
                dimension, SimpleObjPool<AutopoolArray<T> *>{}));
            pool_iter = general_pool_.find(dimension);
        }
        aff_pool_ = &pool_iter->second;
    }

    void require(size_t dimension) {
        if (!aff_pool_) {
            init_aff_pool(dimension);
        }

        dimension_ = dimension;
        if (aff_pool_->cached_list.empty()) {
            data_ = new T[dimension];
            aff_pool_->objects.insert(this);
        } else {
            auto latest_released = aff_pool_->cached_list.top();
            *this = std::move(*latest_released);
            aff_pool_->cached_list.pop();
            // aff_pool_->objects.insert(this);
        }
    }

    void cache() {
        dimension_ = 0;
        if (data_ == nullptr) {
            return;
        }

        aff_pool_->cached_list.push(this);
    }

    void release() {
        dimension_ = 0;
        if (data_ != nullptr) {
            delete[] data_;
        }
    }

private:
    T *data_ = nullptr;

    size_t dimension_ = 0;

    SimpleObjPool<AutopoolArray<T> *> *aff_pool_ = nullptr;
};

template <typename T>
std::map<size_t, SimpleObjPool<AutopoolArray<T> *>>
    AutopoolArray<T>::general_pool_;

} // namespace hehub
