#pragma once

#include "../traits.hpp"

#include <cassert>
#include <optional>
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
  using Node   = typename Ntk::node;
  using signal = typename Ntk::signal;
  using Signal = typename Ntk::signal;

  explicit incomplete_signal_map(Ntk& ntk)
      : ntk_(ntk)
  {
    static_assert(mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type");
    static_assert(traits::has_signal_size_v<Ntk>, "Ntk lacks signal_size()");
    static_assert(traits::has_signal_to_index_v<Ntk>, "Ntk lacks signal_to_index()");

    data_.resize(ntk_.signal_size());
  }

  [[nodiscard]] bool has(signal f) const noexcept {
    auto index = ntk_.signal_to_index(f);
    if (index >= data_.size()) {
      return false;
    }
    return data_[index].has_value();
  }

  void erase(signal f) noexcept {
    auto index = ntk_.signal_to_index(f);
    if (index < data_.size()) {
      data_[index].reset();
    }
  }

  T& operator[](signal f) {
    auto index = ntk_.signal_to_index(f);
    if (index >= data_.size()) {
      data_.resize(index + 1);
    }
    if (!data_[index].has_value()) {
      data_[index].emplace();
    }
    return *data_[index];
  }

  T const& operator[](signal f) const {
    auto index = ntk_.signal_to_index(f);
    assert(index < data_.size() && data_[index].has_value());
    return *data_[index];
  }

  void reset() {
    data_.clear();
    data_.resize(ntk_.signal_size());
  }

  void reset(T const& value) {
    ntk_.foreach_node( [&]( auto const& n ){
      ntk_.foreach_output( n, [&]( auto const& f ){
        data_[f]=value;
      } );
    } );
  }

  void resize() {
    if (ntk_.signal_size() > data_.size()) {
      data_.resize(ntk_.signal_size());
    }
  }

private:
  Ntk& ntk_;
  std::vector<std::optional<T>> data_;
};

} /* namespace network */

} /* namespace rinox */
