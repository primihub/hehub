#pragma once

#include "type_defs.h"
#include <map>
#include <set>
#include <stack>
#include <type_traits>

namespace hehub {

class Poolable {
public:
    Poolable() {}

    virtual void real_release() = 0;
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
        while (!cached_list.empty()) {
            cached_list.pop();
        }

        for (auto obj_ptr : objects) {
            obj_ptr->real_release();
        }
    }
};

template <typename T> class AutopoolArray : public Poolable {
    static std::map<size_t, SimpleObjPool<AutopoolArray<T> *>> general_pool_;

public:
    AutopoolArray() {}

    AutopoolArray(const size_t dimension) {
        require(dimension);
    }

    AutopoolArray(const AutopoolArray &other){
        require(other.dimension_);
        std::copy(other.data_, other.data_ + dimension_, data_);
    }

    AutopoolArray(AutopoolArray &&other) noexcept { *this = std::move(other); }

    ~AutopoolArray() { fake_release(); }

    inline T &operator[](const int idx) { return data_[idx]; }

    inline const T operator[](const int idx) const { return data_[idx]; }

    AutopoolArray &operator=(const AutopoolArray &copying) {
        if (dimension_ != copying.dimension_) {
            fake_release();
            require(copying.dimension_);
        }
        std::copy(copying.data_, copying.data_ + dimension_, data_);
        return *this;
    }

    AutopoolArray &operator=(AutopoolArray &&moving) noexcept {
        if (this == &moving) {
            return *this;
        }

        fake_release();
        dimension_ = moving.dimension_;
        data_ = moving.data_;

        moving.fake_release();
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

    void require(size_t dimension) {
        if (!aff_pool_) {
            auto pool_iter = general_pool_.find(dimension);
            if (pool_iter == general_pool_.end()) {
                general_pool_.insert(std::make_pair(
                    dimension, SimpleObjPool<AutopoolArray<T> *>{}));
                pool_iter = general_pool_.find(dimension);
            }
            aff_pool_ = &pool_iter->second;
        }

        if (aff_pool_->cached_list.empty()) {
            dimension_ = dimension;
            data_ = new T[dimension];
            aff_pool_->objects.insert(this);
        } else {
            auto latest_released = aff_pool_->cached_list.top();
            *this = std::move(*latest_released);
            aff_pool_->cached_list.pop();
            aff_pool_->objects.insert(this);
        }
    }

    void fake_release() {
        dimension_ = 0;
        if (data_ == nullptr) {
            return;
        }

        aff_pool_->cached_list.push(this);
    }

    void real_release() {
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
