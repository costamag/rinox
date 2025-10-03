/* rinox: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file window_manager.hpp
  \brief Construction of windows for mapped networks.

  \author Andrea Costamagna
*/

#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <mockturtle/utils/node_map.hpp>
#include <optional>
#include <queue>
#include <stack>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../math/math.hpp"
#include "../network/visualize.hpp"

namespace rinox {

namespace windowing {

/* -------------------------------------------------------------------------- */
/*                                   color_t                                  */
/* -------------------------------------------------------------------------- */

struct color_t {
  // NOTE: no NONE role — “none” means “not assigned for this window id”
  enum class tag : uint32_t {
    INPS = 1,
    DIVS = 2,
    MFFC = 4,
    TFOS = 8,
    MASK = 0xF
  };

  static constexpr uint32_t num_bits_tag = 4u;  // low bits
  static constexpr uint32_t role_mask = static_cast<uint32_t>(tag::MASK);
  static constexpr uint32_t id_mask = ~role_mask;  // high bits

  uint32_t value{0};

  color_t() = default;
  explicit color_t(tag t) : value(static_cast<uint32_t>(t)) {}

  static color_t from_id(uint32_t id) {
    if (id > (UINT32_MAX >> num_bits_tag))
      throw std::out_of_range("color_t id too large");
    return color_t((id << num_bits_tag), true);
  }

  operator uint32_t() const { return value; }

  /* role manipulation (roles are mutually exclusive) */
  void add(tag t) { value = (value & id_mask) | static_cast<uint32_t>(t); }
  void remove(tag) { value &= id_mask; }  // remove role bits entirely

  bool in(tag t) const {
    uint32_t tag_bits = value & role_mask;
    if (tag_bits == 0) return false;
    if (tag_bits == static_cast<uint32_t>(t)) return true;
    if (tag_bits & static_cast<uint32_t>(t)) {
      throw std::logic_error("Node has incompatible colors");
    }
    return false;
  }

  /* visit id (high bits) */
  uint32_t visited_id() const { return value >> num_bits_tag; }
  void set_visited(uint32_t id) {
    value = (id << num_bits_tag) | (value & role_mask);
  }

 private:
  explicit color_t(uint32_t raw, bool) : value(raw) {}
};

/* -------------------------------------------------------------------------- */
/*                                  window_t                                  */
/* -------------------------------------------------------------------------- */

struct default_window_params_t {
  static constexpr uint32_t inps_capacity = 512u;
  static constexpr uint32_t divs_capacity = 2048u;
  static constexpr uint32_t mffc_capacity = 512u;
  static constexpr uint32_t tfos_capacity = 512u;
  static constexpr uint32_t outs_capacity = 512u;
};

template <class Ntk, typename Params = default_window_params_t>
class window_t {
  struct role_ops {
    std::function<void(typename Ntk::node)> add;
    std::function<void(typename Ntk::node)> remove;
  };

 public:
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;
  using static_params = Params;

  explicit window_t(Ntk& ntk) : ntk_(ntk), map_(ntk) {
    role_map_[color_t::tag::INPS] = {
        [&](node_index_t n) { add_to_inps(n); },
        [&](node_index_t n) { remove_from_inps(n); }};
    role_map_[color_t::tag::DIVS] = {
        [&](node_index_t n) { add_to_divs(n); },
        [&](node_index_t n) { remove_from_divs(n); }};
    role_map_[color_t::tag::MFFC] = {
        [&](node_index_t n) { add_to_mffc(n); },
        [&](node_index_t n) { remove_from_mffc(n); }};
    role_map_[color_t::tag::TFOS] = {
        [&](node_index_t n) { add_to_tfos(n); },
        [&](node_index_t n) { remove_from_tfos(n); }};
  }

  void reset(node_index_t const& n) {
    id_++;
    map_.resize();
    inps_.clear();
    divs_.clear();
    mffc_.clear();
    tfos_.clear();
    outs_.clear();
    pivot = n;
  }

  /* ------------------------------ role changes -----------------------------
   */

  void make(node_index_t const& n, color_t::tag target,
            std::initializer_list<color_t::tag> conflicts) {
    reassign(n, target, conflicts);
  }

  void make_mffc(node_index_t const& n) {
    make(n, color_t::tag::MFFC,
         {color_t::tag::DIVS, color_t::tag::INPS, color_t::tag::TFOS});
  }
  void make_divs(node_index_t const& n) {
    make(n, color_t::tag::DIVS,
         {color_t::tag::MFFC, color_t::tag::INPS, color_t::tag::TFOS});
  }
  void make_inps(node_index_t const& n) {
    make(n, color_t::tag::INPS,
         {color_t::tag::DIVS, color_t::tag::MFFC, color_t::tag::TFOS});
  }
  void make_tfos(node_index_t const& n) {
    make(n, color_t::tag::TFOS,
         {color_t::tag::DIVS, color_t::tag::MFFC, color_t::tag::INPS});
  }

  void make_none(node_index_t const& n) {
    if (map_.has(n) && map_[n].visited_id() == id_) {
      if (map_[n].in(color_t::tag::INPS))
        remove_from_inps(n);
      else if (map_[n].in(color_t::tag::DIVS))
        remove_from_divs(n);
      else if (map_[n].in(color_t::tag::MFFC))
        remove_from_mffc(n);
      else if (map_[n].in(color_t::tag::TFOS))
        remove_from_tfos(n);
      map_[n].remove(color_t::tag::MASK);
    }
  }

  /* --------------------------- marking (seen nodes) ------------------------
   */

  void mark_seen(node_index_t const& n) { map_[n] = color_t::from_id(id_); }
  bool was_seen(node_index_t const& n) const { return is_inside(n); }

  /* ----------------------------- sorting helpers ---------------------------
   */

  void sort_inps() { topological_sort<static_params::inps_capacity>(inps_); }
  void sort_divs() { topological_sort<static_params::divs_capacity>(divs_); }
  void sort_mffc() { topological_sort<static_params::mffc_capacity>(mffc_); }
  void sort_tfos() { topological_sort<static_params::tfos_capacity>(tfos_); }
  void sort_outs() { topological_sort<static_params::outs_capacity>(outs_); }

  /* ------------------------------ membership -------------------------------
   */

  bool is_inps(node_index_t const& n) const {
    return map_.has(n) && (map_[n].visited_id() == id_) &&
           map_[n].in(color_t::tag::INPS);
  }
  bool is_divs(node_index_t const& n) const {
    return map_.has(n) && (map_[n].visited_id() == id_) &&
           map_[n].in(color_t::tag::DIVS);
  }
  bool is_tfos(node_index_t const& n) const {
    return map_.has(n) && (map_[n].visited_id() == id_) &&
           map_[n].in(color_t::tag::TFOS);
  }
  bool is_mffc(node_index_t const& n) const {
    return map_.has(n) && (map_[n].visited_id() == id_) &&
           map_[n].in(color_t::tag::MFFC);
  }
  bool is_none(node_index_t const& n) const {
    return !map_.has(n) || (map_[n].visited_id() != id_);
  }
  bool is_inside(node_index_t const& n) const {
    return map_.has(n) && (map_[n].visited_id() == id_);
  }

  /* -------------------------------- iterate --------------------------------
   */

  template <typename Fn>
  void foreach_input(Fn&& fn) const {
    for (auto i = 0u; i < inps_.size(); i++) fn(inps_[i], i);
  }
  template <typename Fn>
  void foreach_divisor(Fn&& fn) const {
    for (auto i = 0u; i < divs_.size(); i++) fn(divs_[i], i);
  }
  template <typename Fn>
  void foreach_mffc(Fn&& fn) const {
    for (auto i = 0u; i < mffc_.size(); i++) fn(mffc_[i], i);
  }
  template <typename Fn>
  void foreach_tfo(Fn&& fn) const {
    for (auto i = 0u; i < tfos_.size(); i++) fn(tfos_[i], i);
  }
  template <typename Fn>
  void foreach_output(Fn&& fn) const {
    for (auto i = 0u; i < outs_.size(); i++) fn(outs_[i], i);
  }

  /* --------------------------------- debug ---------------------------------
   */

  void print() const {
    std::cout << "INPS" << std::endl;
    for (auto const& l : inps_)
      std::cout << "n[" << ntk_.get_node(l) << "] = ntk.create_pi();\n";
    std::cout << "DIVS" << std::endl;

    for (auto const& d : divs_) {
      if (d.output == 0) {
        auto const ids = ntk_.get_binding_ids(ntk_.get_node(d));
        std::cout << "n[" << ntk_.get_node(d) << "] = ntk.create_node({ ";
        ntk_.foreach_fanin(d, [&](auto const& fi, auto ii) {
          if (ii < (ntk_.fanin_size(d.index) - 1))
            std::cout << "signal{ n[" << fi.index << "], " << fi.output
                      << " }, ";
          else
            std::cout << "signal{ n[" << fi.index << "], " << fi.output
                      << " } }, { ";
        });
        for (auto i = 0u; i + 1u < ids.size(); ++i) std::cout << ids[i] << ", ";
        if (!ids.empty()) std::cout << ids.back();
        std::cout << " }).index;\n";
      }
      std::cout << "ntk.create_po({ n[" << d.index << "], " << d.output
                << " });\n";
    }
    std::cout << "MFFC" << std::endl;

    for (auto const& m : mffc_) {
      auto const ids = ntk_.get_binding_ids(m);
      std::cout << "n[" << m << "] = ntk.create_node({ ";
      ntk_.foreach_fanin(m, [&](auto const& fi, auto ii) {
        if (ii < (ntk_.fanin_size(m) - 1))
          std::cout << "signal{ n[" << fi.index << "], " << fi.output << " }, ";
        else
          std::cout << "signal{ n[" << fi.index << "], " << fi.output
                    << " } }, { ";
      });
      for (auto i = 0u; i + 1u < ids.size(); ++i) std::cout << ids[i] << ", ";
      if (!ids.empty()) std::cout << ids.back();
      std::cout << " }).index;\n";
    }
    std::cout << "TFOS" << std::endl;

    for (auto const& m : tfos_) {
      auto const ids = ntk_.get_binding_ids(m);
      std::cout << "n[" << m << "] = ntk.create_node( ";
      ntk_.foreach_fanin(m, [&](auto const& fi, auto ii) {
        if (ii < (ntk_.fanin_size(m) - 1))
          std::cout << "signal{ n[" << fi.index << "], " << fi.output << " }, ";
        else
          std::cout << "signal{ n[" << fi.index << "], " << fi.output
                    << " } }, { ";
      });
      for (auto i = 0u; i + 1u < ids.size(); ++i) std::cout << ids[i] << ", ";
      if (!ids.empty()) std::cout << ids.back();
      std::cout << " }).index;\n";
    }
    std::cout << "OUTS" << std::endl;

    for (auto const& l : outs_)
      std::cout << "ntk.create_po(signal{ n[" << l.index << "], " << l.output
                << " });\n";
  }

 public:
  node_index_t pivot;
  rinox::math::static_vector<signal_t, static_params::inps_capacity> inps_;
  rinox::math::static_vector<signal_t, static_params::divs_capacity> divs_;
  rinox::math::static_vector<node_index_t, static_params::mffc_capacity> mffc_;
  rinox::math::static_vector<node_index_t, static_params::tfos_capacity> tfos_;
  rinox::math::static_vector<signal_t, static_params::outs_capacity> outs_;

 private:
  /* --------------------------- topological sort ----------------------------
   */

  template <uint32_t Cap>
  void topological_sort(rinox::math::static_vector<node_index_t, Cap>& nodes) {
    topo_++;
    static constexpr uint32_t maskTag = 0xFFFFF000u;    // upper 20 bits
    static constexpr uint32_t maskCount = 0x00000FFFu;  // lower 12 bits
    static constexpr unsigned shiftTag = 12;
    uint32_t tag = topo_;

    // 0) dedup nodes themselves
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());

    // 1) initialize values
    for (auto const& n : nodes) {
      ntk_.set_value(n, (tag << shiftTag) | 0u);
    }

    // 2) compute indegrees (unique fanins only)
    for (auto const& n : nodes) {
      rinox::math::static_vector<node_index_t, Cap> uniq_fanins;
      ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
        auto const ni = ntk_.get_node(fi);
        if ((ntk_.value(ni) >> shiftTag) == tag) {
          if (std::find(uniq_fanins.begin(), uniq_fanins.end(), ni) ==
              uniq_fanins.end())
            uniq_fanins.push_back(ni);
        }
      });
      for (auto const& ni : uniq_fanins) {
        uint32_t cnt = (ntk_.value(n) & maskCount) + 1;
        if (cnt > maskCount)
          throw std::overflow_error(
              "indegree overflow (>4095) in topological_sort");
        ntk_.set_value(n, (tag << shiftTag) | cnt);
      }
    }

    // 3) initialize queue with 0-indegree nodes
    std::queue<node_index_t> q;
    for (auto const& n : nodes) {
      if ((ntk_.value(n) & maskCount) == 0) q.push(n);
    }

    rinox::math::static_vector<node_index_t, Cap> result;

    // 4) Kahn’s algorithm (decrement unique fanouts)
    while (!q.empty()) {
      auto n = q.front();
      q.pop();
      result.push_back(n);

      rinox::math::static_vector<node_index_t, Cap> uniq_fanouts;
      ntk_.foreach_fanout(n, [&](auto const& no) {
        if ((ntk_.value(no) >> shiftTag) == tag) {
          if (std::find(uniq_fanouts.begin(), uniq_fanouts.end(), no) ==
              uniq_fanouts.end())
            uniq_fanouts.push_back(no);
        }
      });

      for (auto const& no : uniq_fanouts) {
        uint32_t cnt = (ntk_.value(no) & maskCount);
        if (cnt == 0)
          throw std::logic_error("indegree underflow in topological_sort");
        cnt--;
        ntk_.set_value(no, (tag << shiftTag) | cnt);
        if (cnt == 0) q.push(no);
      }
    }

    // 5) cycle check
    if (result.size() != nodes.size()) {
      std::cout << "Possible cycle detected in topological_sort\n";
    }

    nodes = std::move(result);
  }

  template <uint32_t Cap>
  void topological_sort(rinox::math::static_vector<signal_t, Cap>& signals) {
    rinox::math::static_vector<node_index_t, Cap> nodes;
    for (auto const& s : signals)
      if (s.output == 0) nodes.push_back(ntk_.get_node(s));

    topological_sort<Cap>(nodes);

    rinox::math::static_vector<signal_t, Cap> result;
    for (auto const& n : nodes)
      ntk_.foreach_output(n, [&](auto const& f) { result.push_back(f); });

    signals = std::move(result);
  }

  void add_to_inps(node_index_t const& n) {
    map_[n] = color_t::from_id(id_);
    map_[n].add(color_t::tag::INPS);
    ntk_.foreach_output(n, [&](auto const& f) {
      divs_.push_back(f);
      inps_.push_back(f);
    });
  }
  void remove_from_inps(node_index_t const& n) {
    if (!map_.has(n)) return;
    map_[n].remove(color_t::tag::INPS);
    ntk_.foreach_output(n, [&](auto const& f) {
      divs_.erase_value(f);
      inps_.erase_value(f);
    });
  }
  void add_to_outs(node_index_t const& n) {
    ntk_.foreach_output(n, [&](auto const& f) { outs_.insert_unique(f); });
  }
  void add_to_divs(node_index_t const& n) {
    map_[n] = color_t::from_id(id_);
    map_[n].add(color_t::tag::DIVS);
    ntk_.foreach_output(n, [&](auto const& f) { divs_.push_back(f); });
  }
  void remove_from_divs(node_index_t const& n) {
    if (!map_.has(n)) return;
    map_[n].remove(color_t::tag::DIVS);
    ntk_.foreach_output(n, [&](auto const& f) { divs_.erase_value(f); });
  }
  void add_to_mffc(node_index_t const& n) {
    map_[n] = color_t::from_id(id_);
    map_[n].add(color_t::tag::MFFC);
    mffc_.push_back(n);
  }
  void remove_from_mffc(node_index_t const& n) {
    if (!map_.has(n)) return;
    map_[n].remove(color_t::tag::MFFC);
    mffc_.erase_value(n);
  }
  void add_to_tfos(node_index_t const& n) {
    map_[n] = color_t::from_id(id_);
    map_[n].add(color_t::tag::TFOS);
    tfos_.push_back(n);
  }
  void remove_from_tfos(node_index_t const& n) {
    if (!map_.has(n)) return;
    map_[n].remove(color_t::tag::TFOS);
    tfos_.erase_value(n);
  }
  void reassign(node_index_t const& n, color_t::tag target,
                std::initializer_list<color_t::tag> conflicts) {
    if (map_.has(n) && map_[n].visited_id() == id_) {
      if (map_[n].in(target)) return;
      for (auto c : conflicts) {
        if (map_[n].in(c)) {
          role_map_[c].remove(n);
        }
      }
    }
    role_map_[target].add(n);
  }

  /* topological_sort methods unchanged (same as you had before) */

 private:
  Ntk& ntk_;
  uint32_t id_ = 0u;
  uint32_t topo_ = 1u;

  mockturtle::incomplete_node_map<color_t, Ntk> map_;
  std::unordered_map<color_t::tag, role_ops> role_map_;
};

/* -------------------------------------------------------------------------- */
/*                            window_manager (public) */
/* -------------------------------------------------------------------------- */

struct default_window_manager_params {
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 10;
  int32_t odc_levels = 0;
  static constexpr uint32_t max_num_leaves = 8;
};

struct window_manager_stats {
  bool valid = false;
};

template <class Ntk, typename Params = default_window_manager_params>
class window_manager {
 public:
  
  using node_index_t = typename Ntk::node;
  using signal_t = typename Ntk::signal;
  static constexpr uint32_t kMaxNumLeaves = Params::max_num_leaves;
  static constexpr uint32_t inps_capacity =
      window_t<Ntk>::static_params::inps_capacity;
  static constexpr uint32_t divs_capacity =
      window_t<Ntk>::static_params::divs_capacity;
  static constexpr uint32_t mffc_capacity =
      window_t<Ntk>::static_params::mffc_capacity;
  static constexpr uint32_t tfos_capacity =
      window_t<Ntk>::static_params::tfos_capacity;
  static constexpr uint32_t outs_capacity =
      window_t<Ntk>::static_params::outs_capacity;

  using input_signals_t = rinox::math::static_vector<signal_t, inps_capacity>;
  using front_signals_t =
      rinox::math::static_vector<signal_t, inps_capacity + outs_capacity>;
  using front_nodes_t = rinox::math::static_vector<node_index_t, divs_capacity>;

 public:
  window_manager(Ntk& ntk, Params const& ps, window_manager_stats& st)
      : ntk_(ntk), window_(ntk_), ps_(ps), st_(st) {}

  bool run(node_index_t const& n) {
    st_.valid = false;

    if (ntk_.fanin_size(n) == 0) return false;
    if (ntk_.is_pi(n)) return false;

    window_.reset(n);

    if (!collect_inps()) return false;
    if (!collect_mffc()) return false;
    if (!collect_divs()) return false;
    if (!collect_tfos()) return false;
    if (!expand_inps()) return false;
    if (!fix_dominos()) return false;

    window_.sort_inps();
    window_.sort_divs();
    window_.sort_mffc();
    window_.sort_tfos();
    window_.sort_outs();

    fix_labels();

    st_.valid = window_.divs_.size() <= ps_.max_num_divisors;
    st_.valid &= window_.inps_.size() <= ps_.max_num_leaves;

    return st_.valid;
  }

  bool is_valid() const { return st_.valid; }

  /* quick queries */
  size_t size() const {
    size_t s = 0;
    s += window_.inps_.size();
    s += window_.divs_.size();
    s += window_.mffc_.size() * Ntk::max_num_outputs;
    s += window_.tfos_.size() * Ntk::max_num_outputs;
    return s;
  }

  bool is_input(node_index_t const& n) const { return window_.is_inps(n); }
  bool is_contained(node_index_t const& n) const {
    return window_.is_inps(n) || window_.is_divs(n) || window_.is_mffc(n) ||
           window_.is_tfos(n);
  }

  size_t num_inputs() const { return window_.inps_.size(); }
  size_t num_outputs() const { return window_.outs_.size(); }
  size_t num_mffc() const { return window_.mffc_.size(); }
  size_t num_divisors() const { return window_.divs_.size(); }

  auto const& get_divisors() const { return window_.divs_; }
  signal_t const& get_divisor(uint32_t const& index) const {
    return window_.divs_[index];
  }
  auto const& get_tfos() const { return window_.tfos_; }
  auto const& get_outputs() const { return window_.outs_; }
  auto const& get_inputs() const { return window_.inps_; }
  auto const& get_divs() const { return window_.divs_; }
  auto const& get_mffc() const { return window_.mffc_; }
  node_index_t const& get_pivot() const { return window_.pivot; }
  window_t<Ntk> const& get_window() const { return window_; }

  template <typename Fn>
  void foreach_input(Fn&& fn) const {
    window_.foreach_input(fn);
  }
  template <typename Fn>
  void foreach_divisor(Fn&& fn) const {
    window_.foreach_divisor(fn);
  }
  template <typename Fn>
  void foreach_mffc(Fn&& fn) const {
    window_.foreach_mffc(fn);
  }
  template <typename Fn>
  void foreach_tfo(Fn&& fn) const {
    window_.foreach_tfo(fn);
  }
  template <typename Fn>
  void foreach_output(Fn&& fn) const {
    window_.foreach_output(fn);
  }

  void print() const { window_.print(); }

  void MarkTfoInside(node_index_t n, uint32_t color)
  {
    if (!is_contained(n)) return;
    ntk_.set_visited(n, color);
    ntk_.foreach_fanout(n, [&](auto const no){
      MarkTfoInside(no, color);
    });
  }


 private:
  /* ----------------------------- INPS collection ---------------------------
   */

  bool collect_inps() {
    front_nodes_t front;
    front_signals_t best_inputs;

    size_t num_inps = 0;
    size_t num_divs = 0;

    window_.mark_seen(window_.pivot);

    ntk_.foreach_fanin(window_.pivot, [&](auto const& fi, auto) {
      auto const ni = ntk_.get_node(fi);
      if (!window_.was_seen(ni)) {
        window_.mark_seen(ni);
        front.push_back(ni);
        num_inps += ntk_.num_outputs(ni);
        num_divs += ntk_.num_outputs(ni);
      }
    });

    best_inputs.clear();
    for (auto const& n : front) {
      ntk_.foreach_output(n, [&](auto const& f) { best_inputs.push_back(f); });
    }

    while (num_inps < (ps_.max_num_leaves + ps_.max_reconvergence)) {
      int best_cost = std::numeric_limits<int>::max();
      std::optional<node_index_t> best_node;

      for (int i = 0; i < static_cast<int>(front.size()); ++i) {
        auto const& n = front[i];
        if (ntk_.is_pi(n)) continue;

        int cost = -static_cast<int>(ntk_.num_outputs(n));
        ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.was_seen(ni)) cost += ntk_.num_outputs(ni);
        });

        if (cost < best_cost) {
          best_cost = cost;
          best_node = n;
        }
      }

      if (best_node &&
          ((num_inps + best_cost) <=
           (ps_.max_num_leaves + ps_.max_reconvergence)) &&
          ((num_divs + best_cost + ntk_.num_outputs(*best_node)) <=
           ps_.max_num_divisors)) {
        ntk_.foreach_fanin(*best_node, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.was_seen(ni)) {
            window_.mark_seen(ni);
            front.push_back(ni);
            num_inps += ntk_.num_outputs(ni);
            num_divs += ntk_.num_outputs(ni);
          }
        });
        num_inps -= ntk_.num_outputs(*best_node);
        front.erase_value(*best_node);

        if (num_inps <= ps_.max_num_leaves) {
          best_inputs.clear();
          for (auto const& n : front) {
            ntk_.foreach_output(
                n, [&](auto const& f) { best_inputs.push_back(f); });
          }
        }
      } else {
        break;
      }
    }

    for (auto const& f : best_inputs) window_.make_inps(ntk_.get_node(f));

    return (window_.inps_.size() <= ps_.max_num_leaves) &&
           (window_.divs_.size() <= ps_.max_num_divisors);
  }

  /* ----------------------------- MFFC collection ---------------------------
   */

  bool collect_mffc() {
    // Root of the cone
    window_.make_mffc(window_.pivot);

    size_t prev_size = 0;
    while (window_.mffc_.size() > prev_size) {
      prev_size = window_.mffc_.size();

      for (size_t i = 0; i < prev_size; ++i) {
        auto const n = window_.mffc_[i];

        ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);

          // Do not absorb the cut boundary or PIs
          if (window_.is_inps(ni)) return;
          if (ntk_.is_pi(ni)) return;
          if (ntk_.is_po(ni)) return;
          // (A fanin cannot be a PO, so no need to check ntk_.is_po(ni))

          // Standard MFFC test: all fanouts already in the current MFFC
          bool all_in = true;
          ntk_.foreach_fanout(
              ni, [&](auto const& fo) { all_in &= window_.is_mffc(fo); });

          if (all_in) {
            window_.make_mffc(ni);
          }
        });
      }

      if (window_.mffc_.size() >= mffc_capacity) break;
    }

    return window_.mffc_.size() <= mffc_capacity;
  }

  /* ------------------------------- expand_inps -----------------------------
   */

  bool expand_inps() {
    while (window_.inps_.size() < ps_.max_num_leaves) {
      int best_cost = std::numeric_limits<int>::max();
      std::optional<node_index_t> best_node;

      const size_t num_inps = window_.inps_.size();
      for (int i = 0; i < static_cast<int>(num_inps); ++i) {
        auto const& f = window_.inps_[i];
        if (f.output != 0) continue;
        auto const n = ntk_.get_node(f);
        if (ntk_.is_pi(n)) continue;

        int cost = -static_cast<int>(ntk_.num_outputs(n));
        ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.is_inps(ni) && !window_.is_divs(ni))
            cost += ntk_.num_outputs(ni);
        });

        if (cost < best_cost) {
          best_cost = cost;
          best_node = n;
        }
      }

      if (!best_node) break;

      if ((window_.inps_.size() + best_cost) <= ps_.max_num_leaves &&
          (window_.divs_.size() + best_cost + ntk_.num_outputs(*best_node)) <=
              ps_.max_num_divisors) {
        window_.make_divs(*best_node);
        ntk_.foreach_fanin(*best_node, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.is_inps(ni) && !window_.is_divs(ni))
            window_.make_inps(ni);
        });
      } else {
        break;
      }
    }

    return (window_.inps_.size() <= ps_.max_num_leaves) &&
           (window_.divs_.size() <= ps_.max_num_divisors);
  }

  /* ------------------------------ divisors (TFIs) --------------------------
   */

  bool collect_tfis() {
    std::stack<node_index_t> stack;

    for (auto const& n : window_.mffc_) {
      ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
        auto const ni = ntk_.get_node(fi);
        if (!window_.was_seen(ni)) return;
        if (window_.is_inps(ni)) return;
        if (window_.is_mffc(ni)) return;
        stack.push(ni);
      });
    }

    while (!stack.empty()) {
      auto n = stack.top();
      stack.pop();

      if (window_.is_divs(n) || window_.is_inps(n) || window_.is_mffc(n))
        continue;

      window_.make_divs(n);

      ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
        auto const ni = ntk_.get_node(fi);
        if (!window_.was_seen(ni)) return;
        if (window_.is_inps(ni)) return;
        if (window_.is_mffc(ni)) return;
        if (window_.is_divs(ni)) return;
        stack.push(ni);
      });
    }

    return window_.divs_.size() <= ps_.max_num_divisors;
  }

  bool collect_side() {
    bool expanded = true;
    while (expanded) {
      expanded = false;

      for (int i = 0; i < static_cast<int>(window_.divs_.size()) &&
                      window_.divs_.size() <= ps_.max_num_divisors;
           ++i) {
        auto const& f = window_.divs_[i];
        if (f.output != 0) continue;

        auto const n = ntk_.get_node(f);

        ntk_.foreach_fanout(n, [&](auto const& no) {
          if (window_.is_mffc(no)) return;
          if (window_.is_divs(no)) return;
          if (window_.is_inps(no)) return;

          bool add = true;
          ntk_.foreach_fanin(no, [&](auto const& fi, auto) {
            auto const ni = ntk_.get_node(fi);
            add &= (window_.is_divs(ni) || window_.is_inps(ni));
          });

          if (add) {
            window_.make_divs(no);
            expanded = true;
          }
        });
      }
    }

    return window_.divs_.size() <= ps_.max_num_divisors;
  }

  bool collect_divs() {
    if (!collect_tfis()) return false;
    if (!collect_side()) return false;
    return (window_.divs_.size() <= ps_.max_num_divisors);
  }

  /* ------------------------------- TFOS growth -----------------------------
   */

  bool collect_tfos() {
    // Case 1: no ODC levels -> pivot is only output
    if (ps_.odc_levels == 0) {
      ntk_.foreach_output(window_.pivot, [&](auto const& f) {
        window_.outs_.insert_unique(f);
      });
      return true;
    }

    // Case 2: we want to explore TFOs
    rinox::math::static_vector<node_index_t, tfos_capacity> current_front;
    rinox::math::static_vector<node_index_t, tfos_capacity> next_front;

    ntk_.incr_trav_id();
    ntk_.set_visited(window_.pivot, ntk_.trav_id());

    // -- add immediate fanouts of pivot
    ntk_.foreach_fanout(window_.pivot, [&](auto const& no) {
      if (ntk_.visited(no) != ntk_.trav_id()) {
        current_front.insert_unique(no);
        ntk_.set_visited(no, ntk_.trav_id());
      }
    });

    // -- count how many NEW inputs are needed
    size_t extra_inputs = 0;
    rinox::math::static_vector<node_index_t, tfos_capacity> extra_nodes;

    for (auto const& n : current_front) {
      ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
        auto const ni = ntk_.get_node(fi);
        if (!window_.is_inside(ni) && ntk_.visited(ni) != ntk_.trav_id()) {
          extra_inputs += ntk_.num_outputs(ni);
          extra_nodes.push_back(ni);
          ntk_.set_visited(ni, ntk_.trav_id());
        }
      });
    }

    // -- check leaf budget
    if ((window_.inps_.size() + extra_inputs) > ps_.max_num_leaves) {
      // Too many inputs: pivot is output instead
      ntk_.foreach_output(window_.pivot, [&](auto const& f) {
        window_.outs_.insert_unique(f);
      });
      return true;
    }

    // -- success: commit the nodes
    for (auto const& n : current_front) window_.make_tfos(n);
    for (auto const& n : extra_nodes) window_.make_inps(n);

    // -- grow up to odc_levels
    for (int level = 1; level < ps_.odc_levels; ++level) {
      next_front.clear();

      for (auto const& n : current_front) {
        ntk_.foreach_fanout(n, [&](auto const& no) {
          if (ntk_.visited(no) != ntk_.trav_id()) {
            next_front.insert_unique(no);
            ntk_.set_visited(no, ntk_.trav_id());
          }
        });
      }

      if (next_front.empty()) break;

      // count new inputs needed
      extra_inputs = 0;
      extra_nodes.clear();

      for (auto const& n : next_front) {
        ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.is_inside(ni) && ntk_.visited(ni) != ntk_.trav_id()) {
            extra_inputs += ntk_.num_outputs(ni);
            extra_nodes.push_back(ni);
            ntk_.set_visited(ni, ntk_.trav_id());
          }
        });
      }

      if ((window_.inps_.size() + extra_inputs) > ps_.max_num_leaves) {
        // stop growing here
        break;
      }

      // commit this level
      for (auto const& n : next_front) window_.make_tfos(n);
      for (auto const& n : extra_nodes) window_.make_inps(n);

      current_front = next_front;
    }

    if (window_.tfos_.size() > 0) {
      // -- mark outputs: any TFO that leaves window or is a PO
      window_.foreach_tfo([&](auto const& n, auto) {
        bool output = ntk_.is_po(n);
        ntk_.foreach_fanout(n, [&](auto const& no) {
          if (!window_.is_inside(no)) output = true;
        });
        if (output) {
          ntk_.foreach_output(
              n, [&](auto const& f) { window_.outs_.insert_unique(f); });
        }
      });
    } else {
      ntk_.foreach_output(window_.pivot, [&](auto const& f) {
        window_.outs_.insert_unique(f);
      });
    }

    return (window_.inps_.size() <= ps_.max_num_leaves) &&
           (window_.divs_.size() <= ps_.max_num_divisors) &&
           (window_.tfos_.size() <= tfos_capacity) &&
           (window_.outs_.size() <= outs_capacity);
  }

  bool fix_dominos() {
    bool any_change = false;
    bool changed = true;

    while (changed) {
      changed = false;

      /* Step 1: refine INPS */
      for (int i = 0; i < static_cast<int>(window_.inps_.size()); ++i) {
        auto const f = window_.inps_[i];
        if (f.output != 0) continue;
        auto const n = ntk_.get_node(f);
        if (ntk_.is_pi(n)) continue;

        // Check if all fanins are inside
        bool all_fanins_inside = true;
        ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.is_inps(ni) && !window_.is_divs(ni))
            all_fanins_inside = false;
        });
        if (!all_fanins_inside) continue;

        // Check if all fanouts are MFFC
        window_.make_divs(n);
        changed = true;
        any_change = true;
        break;  // restart since containers have changed
      }
    }

    return (window_.inps_.size() <= ps_.max_num_leaves) &&
           (window_.divs_.size() <= ps_.max_num_divisors) &&
           (window_.tfos_.size() <= tfos_capacity) &&
           (window_.outs_.size() <= outs_capacity);
  }

  void fix_labels()
  {
    std::set<signal_t> inps;
    foreach_input( [&]( auto const& f, auto i ){
      window_.make_inps(ntk_.get_node(f));
      inps.insert(f);
    } );
    foreach_divisor( [&]( auto const& f, auto i ){
      if ( inps.find(f) == inps.end() )
        window_.make_divs(ntk_.get_node(f));
    } );
  }

 private:
  Ntk& ntk_;
  window_t<Ntk> window_;
  Params const& ps_;
  window_manager_stats& st_;
};

}  // namespace windowing
}  // namespace rinox
