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
  \file signal_map.hpp
  \brief Map indexed by network signals

  \author Andrea Costamagna
*/

#pragma once

#include "../traits.hpp"

#include <cassert>
#include <memory>
#include <unordered_map>
#include <variant>
#include <vector>

#include <mockturtle/traits.hpp>

namespace rinox
{

namespace network
{

/*! \brief Vector-based signal map with validity query
 *
 * This container is a variant of the `incomplete_node_map` tailored for cases
 * in which different signals pointing to the same node should store different
 * values. A crucial use case is when signals contain bitfields for specifying
 * the output in a multiple output gate.
 *
 * The implementation uses a vector as underlying data structure, so
 * that it benefits from fast access. It is supplemented with an
 * additional validity field such that it can be used like an
 * `unordered_signal_map`.
 *
 * **Required network functions:**
 * - `signal_size`
 * - `signal_to_index`
 *
 */

template<class T, class Ntk>
class incomplete_signal_map {
public:
  using node   = typename Ntk::node;
  using signal = typename Ntk::signal;
  using container_type = std::vector<std::variant<std::monostate, T>>;

  explicit incomplete_signal_map(Ntk const& ntk) noexcept
      : ntk_(&ntk), data_(std::make_shared<container_type>(ntk.signal_size())) {
    static_assert(!std::is_same_v<T, std::monostate>, "T cannot be std::monostate");
    static_assert(mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type");
    static_assert(traits::has_signal_size_v<Ntk>, "Ntk lacks signal_size()");
    static_assert(traits::has_signal_to_index_v<Ntk>, "Ntk lacks signal_to_index()");
  }

  incomplete_signal_map(Ntk const& ntk, T const& init_value)
      : ntk_(&ntk),
        data_(std::make_shared<container_type>(ntk.signal_size(), init_value)) {
    static_assert(!std::is_same_v<T, std::monostate>, "T cannot be std::monostate");
    static_assert(mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type");
    static_assert(traits::has_signal_size_v<Ntk>, "Ntk lacks signal_size()");
    static_assert(traits::has_signal_to_index_v<Ntk>, "Ntk lacks signal_to_index()");
  }

  [[nodiscard]] auto size() const noexcept -> std::size_t {
    return data_->size();
  }

  [[nodiscard]] bool has(signal const& f) const noexcept {
    assert(ntk_ && "incomplete_signal_map constructed with null network pointer");
    const auto idx = ntk_->signal_to_index(f);
    if (idx >= data_->size()) return false;
    return std::holds_alternative<T>((*data_)[idx]);
  }

  void erase(signal const& f) noexcept {
    const auto idx = ntk_->signal_to_index(f);
    if (idx >= data_->size()) return;
    (*data_)[idx] = std::monostate();
  }

  T& operator[](signal const& f) {
    const auto idx = ntk_->signal_to_index(f);
    if (idx >= data_->size()) data_->resize(idx + 1);
    if (!std::holds_alternative<T>((*data_)[idx])) {
      (*data_)[idx] = T{};
    }
    return std::get<T>((*data_)[idx]);
  }

  T const& operator[](signal const& f) const {
    const auto idx = ntk_->signal_to_index(f);
    assert(idx < data_->size() && "index out of bounds");
    assert(std::holds_alternative<T>((*data_)[idx]) && "uninitialized value");
    return std::get<T>((*data_)[idx]);
  }

  void reset() {
    data_->clear();
    data_->resize(ntk_->signal_size());
  }

  void reset(T const& init_value) {
    data_->clear();
    data_->resize(ntk_->signal_size(), init_value);
  }

  void resize() {
    if (ntk_->signal_size() > data_->size()) {
      data_->resize(ntk_->signal_size());
    }
  }

private:
  Ntk const* ntk_;
  std::shared_ptr<container_type> data_;
};



} /* namespace network */

} /* namespace rinox */