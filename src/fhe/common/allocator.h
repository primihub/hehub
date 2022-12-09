#pragma once

#include "type_defs.h"
#include <cassert>
#include <map>
#include <set>
#include <stack>
#include <type_traits>

namespace hehub {

template <typename T> class FixedBlockAllocator {
public:
    /// Constructor
    /// @param[in] size - size of the fixed blocks
    FixedBlockAllocator(size_t size) : block_size_(size) {}

    /// Destructor
    ~FixedBlockAllocator() {
        while (head_)
            delete[](T *) _pop();
    }

    /// Get a pointer to a memory block.
    /// @return Returns pointer to the block. Otherwise nullptr if
    /// unsuccessful.
    void *allocate() {
        // If can't obtain existing block then get a new one
        void *block = _pop();
        if (!block) {
            blocks_total_++;
            block = (void *)new T[block_size_];
        } else {
            blocks_free_--;
        }

        blocks_in_use_++;

        return block;
    }

    /// Return a pointer to the memory pool.
    /// @param[in] to_cache - block of memory deallocate (i.e push onto
    /// free-list)
    void deallocate(void *to_cache) {
        _push(to_cache);
        blocks_in_use_--;
        blocks_free_++;
    }

    /// Gets the fixed block memory size, in bytes, handled by the allocator.
    /// @return The fixed block size in bytes.
    const size_t get_block_size() const { return block_size_; }

    /// Gets the number of blocks in use.
    /// @return The number of blocks in use by the application.
    const size_t get_blocks_in_use() const { return blocks_in_use_; }

    /// Gets the total number of allocations for this allocator instance.
    /// @return The total number of allocations.
    const size_t get_blocks_total() const { return blocks_total_; }

    /// Gets the total number of deallocations for this allocator instance.
    /// @return The total number of deallocations.
    const size_t get_blocks_free() const { return blocks_free_; }

private:
    /// Push a memory block onto head of free-list.
    /// @param[in] to_cache - block of memory to push onto free-list
    void _push(void *to_cache) {
        Block *block = (Block *)to_cache;
        block->next_ = head_;
        head_ = block;
    }

    /// Pop a memory block from head of free-list.
    /// @return Returns pointer to the block. Otherwise NULL if
    /// unsuccessful.
    void *_pop() {
        Block *block = NULL;

        if (head_) {
            block = head_;
            head_ = head_->next_;
        }

        return (void *)block;
    }

    struct Block {
        Block *next_;
    };

    const size_t block_size_;

    Block *head_ = nullptr;

    size_t blocks_in_use_ = 0;

    size_t blocks_total_ = 0;

    size_t blocks_free_ = 0;
};

template <typename T> class SmartArray {
public:
    SmartArray() {}

    SmartArray(const size_t dimension) : dimension_(dimension) {
        require(dimension_);
    }

    SmartArray(const SmartArray &other) {
        dimension_ = other.dimension_;
        aff_allocator_ = other.aff_allocator_;
        require(dimension_);
        std::copy(other.data_, other.data_ + dimension_, data_);
    }

    SmartArray(SmartArray &&other) noexcept { *this = std::move(other); }

    ~SmartArray() { cache(); }

    inline T &operator[](const int idx) { return data_[idx]; }

    inline const T operator[](const int idx) const { return data_[idx]; }

    SmartArray &operator=(const SmartArray &copying) {
        if (dimension_ != copying.dimension_) {
            cache();
            require(copying.dimension_);
        }
        std::copy(copying.data_, copying.data_ + dimension_, data_);
        return *this;
    }

    SmartArray &operator=(SmartArray &&moving) noexcept {
        if (this == &moving) {
            return *this;
        }

        dimension_ = moving.dimension_;
        moving.dimension_ = 0;

        aff_allocator_ = moving.aff_allocator_;
        if (!aff_allocator_) {
            init_allocator(dimension_);
        }
        moving.aff_allocator_ = nullptr;

        data_ = moving.data_;
        moving.data_ = nullptr;

        return *this;
    }

    inline bool operator==(const SmartArray &other) const {
        if (dimension_ != other.dimension_)
            return false;
        for (int i = 0; i < dimension_; i++) {
            if ((*this)[i] != other[i]) {
                return false;
            }
        }
        return true;
    }

    inline bool operator!=(const SmartArray &comparing) const {
        return !((*this) == comparing);
    }

    inline T *data() { return data_; }

    inline const T *data() const { return data_; }

    inline T *begin() { return data_; }

    inline const T *begin() const { return data_; }

    inline T *end() { return data_ + dimension_; }

    inline const T *end() const { return data_ + dimension_; }

    inline const auto &aff_allocator() const { return *aff_allocator_; }

    void init_allocator(size_t dimension) {
        auto allocator_iter = allocator_hub_.find(dimension);
        if (allocator_iter == allocator_hub_.end()) {
            allocator_hub_.insert(
                std::make_pair(dimension, FixedBlockAllocator<T>(dimension)));
            allocator_iter = allocator_hub_.find(dimension);
        }
        aff_allocator_ = &allocator_iter->second;
    }

    inline void require(size_t dimension) {
        if (!aff_allocator_) {
            init_allocator(dimension_);
        }

        data_ = (T *)aff_allocator_->allocate();
    }

    inline void cache() {
        dimension_ = 0;
        if (aff_allocator_) {
            aff_allocator_->deallocate((void *)data_);
            aff_allocator_ = nullptr;
        }
    }

private:
    static std::map<size_t, FixedBlockAllocator<T>> allocator_hub_;

    T *data_ = nullptr;

    size_t dimension_ = 0;

    FixedBlockAllocator<T> *aff_allocator_ = nullptr;
};

template <typename T>
std::map<size_t, FixedBlockAllocator<T>> SmartArray<T>::allocator_hub_;

} // namespace hehub
