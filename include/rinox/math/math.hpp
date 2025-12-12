/* rinox: C++ logic network library
 * Copyright (C) 2025 EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file constexpr_functions.hpp
  \brief Various functions for compile time evaluation

  \author Andrea Costamagna
*/

#pragma once
#include <assert.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace rinox {

namespace math {

constexpr uint32_t log2_ceil(uint32_t x) {
  uint32_t r = 0;
  uint32_t v = x - 1;
  while (v > 0) {
    ++r;
    v >>= 1;
  }
  return r;
}

template <typename T, std::size_t Cap>
struct static_vector {
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using iterator = T*;
  using const_iterator = const T*;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  static_vector() = default;

  static_vector(std::vector<T> const& vec) {
    for (auto const& x : vec) push_back(x);
  }

  explicit static_vector(size_t count, T const& value = T{}) {
    while (n < count) push_back(value);
  }

  template <class InputIt>
  static_vector(InputIt first, InputIt last) {
    for (; first != last; ++first) push_back(*first);
  }

  static_vector(std::initializer_list<T> init) {
    for (auto const& x : init) push_back(x);
  }

  static_vector(static_vector const& other) {
    for (auto const& x : other) push_back(x);
  }

  static_vector(static_vector&& other) noexcept {
    for (auto& x : other) push_back(std::move(x));
    other.clear();
  }

  static_vector& operator=(static_vector const& other) {
    if (this != &other) {
      clear();
      for (auto const& x : other) push_back(x);
    }
    return *this;
  }

  static_vector& operator=(static_vector&& other) noexcept {
    if (this != &other) {
      clear();
      for (auto& x : other) push_back(std::move(x));
      other.clear();
    }
    return *this;
  }

  ~static_vector() { clear(); }

  std::array<T, Cap> data_;
  std::size_t n = 0;

  constexpr std::size_t size() const noexcept { return n; }
  constexpr std::size_t capacity() const noexcept { return Cap; }
  constexpr bool empty() const noexcept { return n == 0; }
  void clear() noexcept { n = 0; }

  T& operator[](std::size_t i) noexcept { return data_[i]; }
  const T& operator[](std::size_t i) const noexcept { return data_[i]; }
  T& back() {
    if (n == 0) throw std::out_of_range("back on empty static_vector");
    return data_[n - 1];
  }
  const T& back() const {
    if (n == 0) throw std::out_of_range("back on empty static_vector");
    return data_[n - 1];
  }

  void push_back(T const& x) {
    if (n >= Cap) {
      std::cout << n << std::endl;
      assert(n < Cap && "static_vector full");
    }
    data_[n++] = x;
  }

  void push_back(T&& x) {
    if (n >= Cap) {
      assert(false && "!");
      throw std::out_of_range("static_vector full");
    }
    data_[n++] = std::move(x);
  }

  template <typename... Args>
  T& emplace_back(Args&&... args) {
    if (n >= Cap) {
      std::cout << n << std::endl;
      throw std::out_of_range("static_vector full");
    }
    return data_[n++] = T(std::forward<Args>(args)...);
  }

  void pop_back() {
    if (n == 0) throw std::out_of_range("pop_back on empty static_vector");
    --n;
  }

  void erase(std::size_t i) {
    if (i >= n) throw std::out_of_range("erase out of range");
    if (i < n - 1) std::swap(data_[i], data_[n - 1]);
    --n;
  }

  iterator erase(iterator pos) {
    if (pos < begin() || pos >= end())
      throw std::out_of_range("erase iterator out of range");

    size_t idx = static_cast<size_t>(pos - begin());

    // shift left
    for (size_t i = idx; i + 1 < n; ++i) data_[i] = std::move(data_[i + 1]);

    --n;
    return begin() + idx;  // return iterator to next element
  }

  iterator erase(iterator first, iterator last) {
    if (first < begin() || last > end() || first > last)
      throw std::out_of_range("erase range out of range");

    size_t idx_first = static_cast<size_t>(first - begin());
    size_t idx_last = static_cast<size_t>(last - begin());
    size_t count = idx_last - idx_first;

    for (size_t i = idx_first; i + count < n; ++i)
      data_[i] = std::move(data_[i + count]);

    n -= count;
    return begin() + idx_first;  // iterator to next valid element
  }

  void erase_value(T const& v) {
    if (empty()) return;
    for (int i = n - 1; i >= 0; i--)
      if (data_[i] == v) erase(i);
  }

  void resize(size_t new_size, T const& v) {
    if (new_size > Cap) {
      assert(false && "resize");
      throw std::out_of_range("resize out of range");
    }
    while (n < new_size) push_back(v);
  }

  void insert_unique(T const& v) {
    if (n >= Cap) {
      std::cout << n << std::endl;
      throw std::out_of_range("vector full");
    }
    for (std::size_t i = 0; i < n; ++i)
      if (data_[i] == v) return;
    push_back(v);
  }

  iterator insert(iterator pos, T const& value) {
    if (n >= Cap) {
      std::cout << n << std::endl;
      throw std::out_of_range("static_vector full");
    }

    // compute index from iterator
    size_t idx = static_cast<size_t>(pos - begin());

    // shift elements one step to the right
    for (size_t i = n; i > idx; --i) data_[i] = std::move(data_[i - 1]);

    data_[idx] = value;
    ++n;
    return begin() + idx;
  }

  iterator insert(iterator pos, T&& value) {
    if (n >= Cap) throw std::out_of_range("static_vector full");

    size_t idx = static_cast<size_t>(pos - begin());

    for (size_t i = n; i > idx; --i) data_[i] = std::move(data_[i - 1]);

    data_[idx] = std::move(value);
    ++n;
    return begin() + idx;
  }

  void resize(size_t new_size) {
    if (new_size > Cap) {
      assert(false && "resize");
      throw std::out_of_range("resize out of range");
    }
    n = new_size;
  }

  void shift_insert(size_t j, T const& v) {
    if (n >= Cap) throw std::out_of_range("shif insert out of range");
    for (auto i = n; i > j; --i) data_[i] = data_[i - 1];
    data_[j] = v;
    n++;
  }

  bool operator==(static_vector<T, Cap> const& other) const {
    if (n != other.n) return false;
    for (int i = 0; i < n; ++i) {
      if (data_[i] != other.data_[i]) return false;
    }
    return true;
  }

  bool operator==(std::vector<T> const& other) const {
    if (size() != other.size()) return false;
    for (int i = 0; i < n; ++i)
      if (data_[i] != other[i]) return false;
    return true;
  }

  bool operator!=(std::vector<T> const& other) const {
    if (size() != other.size()) return true;
    for (int i = 0; i < n; ++i)
      if (data_[i] != other[i]) return true;
    return false;
  }

  bool operator!=(static_vector<T, Cap> const& other) const {
    return (n != other.n) || (data_ != other.data_);
  }

  T* begin() noexcept { return data_.data(); }
  T* end() noexcept { return data_.data() + n; }
  const T* begin() const noexcept { return data_.data(); }
  const T* end() const noexcept { return data_.data() + n; }
  const T* const cbegin() const noexcept { return data_.data(); }
  const T* const cend() const noexcept { return data_.data() + n; }
};

}  // namespace math

} /* namespace rinox */