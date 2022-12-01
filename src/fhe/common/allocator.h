#pragma once

#include "type_defs.h"
#include <cassert>
#include <map>
#include <set>
#include <stack>
#include <type_traits>

namespace hehub {

// class Poolable {
// public:
//     Poolable() {}

//     virtual void release() = 0;
// };

// template <typename T> class AutopoolArray : public Poolable {
template <typename T> class AutopoolArray {
public:
    struct Allocator {
        std::stack<T *> cached_list;

        std::set<T *> arrays;

        ~Allocator() {
            std::stack<T *>().swap(cached_list);

            for (auto obj_ptr : arrays) {
                delete[] obj_ptr;
            }
            std::set<T *>().swap(arrays);
        }
    };
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

        allocator_ = moving.allocator_;
        if (!allocator_) {
            init_aff_pool(moving.dimension_);
        }
        if (allocator_->arrays.find(data_) != allocator_->arrays.end()) {
            cache();
        }

        assert(allocator_->arrays.find(moving.data_) !=
               allocator_->arrays.end());
        dimension_ = moving.dimension_;
        moving.dimension_ = 0;
        data_ = moving.data_;
        moving.data_ = nullptr;

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

    inline const auto &aff_pool() const { return *allocator_; }

    void init_aff_pool(size_t dimension) {
        auto pool_iter = general_pool_.find(dimension);
        if (pool_iter == general_pool_.end()) {
            general_pool_.insert(std::make_pair(dimension, Allocator{}));
            pool_iter = general_pool_.find(dimension);
        }
        allocator_ = &pool_iter->second;
    }

    void require(size_t dimension) {
        if (!allocator_) {
            init_aff_pool(dimension);
        }

        dimension_ = dimension;
        if (allocator_->cached_list.empty()) {
            data_ = new T[dimension];
            allocator_->arrays.insert(data_);
        } else {
            auto latest_released = allocator_->cached_list.top();
            data_ = latest_released;
            allocator_->cached_list.pop();
        }
    }

    void cache() {
        dimension_ = 0;
        if (data_ == nullptr) {
            allocator_ = nullptr;
            return;
        }

        allocator_->cached_list.push(data_);
        allocator_ = nullptr;
    }

private:
    static std::map<size_t, Allocator> general_pool_;

    T *data_ = nullptr;

    size_t dimension_ = 0;

    Allocator *allocator_ = nullptr;
};

template <typename T>
std::map<size_t, typename AutopoolArray<T>::Allocator>
    AutopoolArray<T>::general_pool_;

} // namespace hehub
