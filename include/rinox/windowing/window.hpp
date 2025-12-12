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
  \file struct.hpp
  \brief Structural resinthesis of LUT networks

  \author Andrea Costamagna
*/

#ifndef RINOX_WINDOWING_WINDOW_HPP_
#define RINOX_WINDOWING_WINDOW_HPP_

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

#include "../analyzers/evaluators/times_evaluator.hpp"
#include "../math/math.hpp"
#include "../network/visualize.hpp"
#include "window_tags.hpp"

namespace rinox {

namespace windowing {

struct WindowParamsDefault {
  static constexpr uint32_t kInpsCapacity = 512u;
  static constexpr uint32_t kDivsCapacity = 2048u;
  static constexpr uint32_t kMffcCapacity = 512u;
  static constexpr uint32_t kTfosCapacity = 512u;
  static constexpr uint32_t kOutsCapacity = 512u;
};

struct WindowStatsDefault {};

template <class Ntk, typename Params = WindowParamsDefault,
          typename Stats = WindowStatsDefault>
class Window {
  struct RoleOps {
    std::function<void(typename Ntk::node)> Add;
    std::function<void(typename Ntk::node)> Remove;
  };

 public:
  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using StaticParams = Params;

  explicit Window(Ntk& ntk) : ntk_(ntk), map_(ntk) {
    role_map_[WindowTags::tag::INPS] = {[&](Node n) { AddToInps(n); },
                                        [&](Node n) { RemoveFromInps(n); }};
    role_map_[WindowTags::tag::DIVS] = {[&](Node n) { AddToDivs(n); },
                                        [&](Node n) { RemoveFromDivs(n); }};
    role_map_[WindowTags::tag::MFFC] = {[&](Node n) { AddToMffc(n); },
                                        [&](Node n) { RemoveFromMffc(n); }};
    role_map_[WindowTags::tag::TFOS] = {[&](Node n) { AddToTfos(n); },
                                        [&](Node n) { RemoveFromTfos(n); }};
  }

  void Reset(Node const& n) {
    id_++;
    map_.resize();
    inps_.clear();
    divs_.clear();
    mffc_.clear();
    tfos_.clear();
    outs_.clear();
    root_ = n;
  }

  void Make(Node const& n, WindowTags::tag target,
            std::initializer_list<WindowTags::tag> conflicts) {
    Reassign(n, target, conflicts);
  }

  void MakeMffc(Node const& n) {
    Make(n, WindowTags::tag::MFFC,
         {WindowTags::tag::DIVS, WindowTags::tag::INPS, WindowTags::tag::TFOS});
  }
  void MakeDivs(Node const& n) {
    Make(n, WindowTags::tag::DIVS,
         {WindowTags::tag::MFFC, WindowTags::tag::INPS, WindowTags::tag::TFOS});
  }
  void MakeInps(Node const& n) {
    Make(n, WindowTags::tag::INPS,
         {WindowTags::tag::DIVS, WindowTags::tag::MFFC, WindowTags::tag::TFOS});
  }
  void MakeTfos(Node const& n) {
    Make(n, WindowTags::tag::TFOS,
         {WindowTags::tag::DIVS, WindowTags::tag::MFFC, WindowTags::tag::INPS});
  }

  void MakeNone(Node const& n) {
    if (map_.has(n) && map_[n].VisitedId() == id_) {
      if (map_[n].In(WindowTags::tag::INPS))
        RemoveFromInps(n);
      else if (map_[n].In(WindowTags::tag::DIVS))
        RemoveFromDivs(n);
      else if (map_[n].In(WindowTags::tag::MFFC))
        RemoveFromMffc(n);
      else if (map_[n].In(WindowTags::tag::TFOS))
        RemoveFromTfos(n);
      map_[n].Remove(WindowTags::tag::MASK);
    }
  }

  /* --------------------------------------------------------------------- */
  /*                           counting helpers                            */
  /* --------------------------------------------------------------------- */
  size_t NumInps() const { return inps_.size(); }
  size_t NumTfos() const { return tfos_.size(); }
  size_t NumOuts() const { return outs_.size(); }
  size_t NumMffc() const { return mffc_.size(); }
  size_t NumDivs() const { return divs_.size(); }

  /* --------------------------------------------------------------------- */
  /*                          marking seen helpers                         */
  /* --------------------------------------------------------------------- */
  void MarkSeen(Node const& n) { map_[n] = WindowTags::FromId(id_); }
  bool WasSeen(Node const& n) const { return IsInside(n); }

  /* --------------------------------------------------------------------- */
  /*                           ownership helpers                           */
  /* --------------------------------------------------------------------- */
  bool IsInps(Node const& n) const {
    return map_.has(n) && (map_[n].VisitedId() == id_) &&
           map_[n].In(WindowTags::tag::INPS);
  }
  bool IsDivs(Node const& n) const {
    return map_.has(n) && (map_[n].VisitedId() == id_) &&
           map_[n].In(WindowTags::tag::DIVS);
  }
  bool IsTfos(Node const& n) const {
    return map_.has(n) && (map_[n].VisitedId() == id_) &&
           map_[n].In(WindowTags::tag::TFOS);
  }
  bool IsMffc(Node const& n) const {
    return map_.has(n) && (map_[n].VisitedId() == id_) &&
           map_[n].In(WindowTags::tag::MFFC);
  }
  bool IsNone(Node const& n) const {
    return !map_.has(n) || (map_[n].VisitedId() != id_);
  }
  bool IsInside(Node const& n) const {
    return map_.has(n) && (map_[n].VisitedId() == id_);
  }

  /* --------------------------------------------------------------------- */
  /*                                  iterators                            */
  /* --------------------------------------------------------------------- */

  template <typename Fn>
  void ForeachNodeMffc(Fn&& fn) const {
    for (auto i = 0u; i < mffc_.size(); i++) fn(mffc_[i], i);
  }
  template <typename Fn>
  void ForeachNodeTfos(Fn&& fn) const {
    for (auto i = 0u; i < tfos_.size(); i++) fn(tfos_[i], i);
  }
  template <typename Fn>
  void ForeachSignalOuts(Fn&& fn) const {
    for (auto i = 0u; i < outs_.size(); i++) fn(outs_[i], i);
  }
  template <typename Fn>
  void ForeachSignalInps(Fn&& fn) const {
    for (auto i = 0u; i < inps_.size(); i++) fn(inps_[i], i);
  }
  template <typename Fn>
  void ForeachSignalDivs(Fn&& fn) const {
    for (auto i = 0u; i < divs_.size(); i++) fn(divs_[i], i);
  }

  void AddToInps(Node const& n) {
    map_[n] = WindowTags::FromId(id_);
    map_[n].Add(WindowTags::tag::INPS);
    ntk_.foreach_output(n, [&](auto const& f) {
      divs_.push_back(f);
      inps_.push_back(f);
    });
  }
  void RemoveFromInps(Node const& n) {
    if (!map_.has(n)) return;
    map_[n].Remove(WindowTags::tag::INPS);
    ntk_.foreach_output(n, [&](auto const& f) {
      divs_.erase_value(f);
      inps_.erase_value(f);
    });
  }
  void AddToOuts(Node const& n) {
    ntk_.foreach_output(n, [&](auto const& f) { outs_.insert_unique(f); });
  }
  void AddToDivs(Node const& n) {
    map_[n] = WindowTags::FromId(id_);
    map_[n].Add(WindowTags::tag::DIVS);
    ntk_.foreach_output(n, [&](auto const& f) { divs_.push_back(f); });
  }
  void RemoveFromDivs(Node const& n) {
    if (!map_.has(n)) return;
    map_[n].Remove(WindowTags::tag::DIVS);
    ntk_.foreach_output(n, [&](auto const& f) { divs_.erase_value(f); });
  }
  void AddToMffc(Node const& n) {
    map_[n] = WindowTags::FromId(id_);
    map_[n].Add(WindowTags::tag::MFFC);
    mffc_.push_back(n);
  }
  void RemoveFromMffc(Node const& n) {
    if (!map_.has(n)) return;
    map_[n].Remove(WindowTags::tag::MFFC);
    mffc_.erase_value(n);
  }
  void AddToTfos(Node const& n) {
    map_[n] = WindowTags::FromId(id_);
    map_[n].Add(WindowTags::tag::TFOS);
    tfos_.push_back(n);
  }
  void RemoveFromTfos(Node const& n) {
    if (!map_.has(n)) return;
    map_[n].Remove(WindowTags::tag::TFOS);
    tfos_.erase_value(n);
  }
  void Reassign(Node const& n, WindowTags::tag target,
                std::initializer_list<WindowTags::tag> conflicts) {
    if (map_.has(n) && map_[n].VisitedId() == id_) {
      if (map_[n].In(target)) return;
      for (auto c : conflicts) {
        if (map_[n].In(c)) {
          role_map_[c].Remove(n);
        }
      }
    }
    role_map_[target].Add(n);
  }

 public:
  auto const& GetSignalsInps() { return inps_; }
  auto const& GetSignalsOuts() { return outs_; }
  auto const& GetSignalsDivs() { return divs_; }
  auto const& GetNodesMffc() { return mffc_; }
  auto const& GetNodesTfos() { return tfos_; }
  Node GetRoot() { return root_; }

 public:
  Node root_;
  rinox::math::static_vector<Signal, Params::kInpsCapacity> inps_;
  rinox::math::static_vector<Signal, Params::kDivsCapacity> divs_;
  rinox::math::static_vector<Node, Params::kMffcCapacity> mffc_;
  rinox::math::static_vector<Node, Params::kTfosCapacity> tfos_;
  rinox::math::static_vector<Signal, Params::kOutsCapacity> outs_;

 private:
  Ntk& ntk_;
  uint32_t id_ = 0u;
  mockturtle::incomplete_node_map<WindowTags, Ntk> map_;
  std::unordered_map<WindowTags::tag, RoleOps> role_map_;
};

}  // namespace windowing
}  // namespace rinox

#endif