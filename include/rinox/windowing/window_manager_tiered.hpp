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

// Structural resynthesis window manager for LUT networks.

#ifndef RINOX_WINDOWING_WINDOW_MANAGER_HPP_
#define RINOX_WINDOWING_WINDOW_MANAGER_HPP_

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
#include "window.hpp"

namespace rinox {
namespace windowing {

// Default parameters for WindowManager. Controls window size and depth limits.
// Inherit from the raw window parameters to share the common information.
struct WindowManagerParamsDefault : WindowParamsDefault {
  static constexpr uint32_t kInpsCapacity = 512u;
  static constexpr uint32_t kDivsCapacity = 2048u;
  static constexpr uint32_t kMffcCapacity = 512u;
  static constexpr uint32_t kTfosCapacity = 512u;
  static constexpr uint32_t kOutsCapacity = 512u;
  static constexpr uint32_t kMaxNumLeaves = 8u;
  uint32_t max_num_divisors = 128u;
  uint32_t max_reconvergence = 10u;
  int max_num_levels_tfi = 10;
  int max_num_levels_tfo = 0;
  bool collect_sides = true;

  void Report() {
    std::cout << "Window Manager Parameters:" << std::endl;
    std::cout << "  - Max number of leaves: " << kMaxNumLeaves << std::endl;
    std::cout << "  - Max number of divisors: " << max_num_divisors
              << std::endl;
    std::cout << "  - Max reconvergence: " << max_reconvergence << std::endl;
    std::cout << "  - Max TFI levels: " << max_num_levels_tfi << std::endl;
    std::cout << "  - Max TFO levels: " << max_num_levels_tfo << std::endl;
  }
};

// Default statistics structure for WindowManager.
struct WindowManagerStatsDefault : WindowStatsDefault {
  bool valid = false;
};

// Manages extraction of logic windows around nodes in a network.
template <class Ntk, typename Params = WindowManagerParamsDefault,
          typename Stats = WindowManagerStatsDefault>
class WindowManager {
 public:
  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using WindowT = Window<Ntk, Params, Stats>;
  using TimesEvaluator = typename analyzers::evaluators::TimesEvaluator<Ntk>;

  using FrontSignals =
      rinox::math::static_vector<Signal,
                                 Params::kInpsCapacity + Params::kOutsCapacity>;
  using FrontNodes = rinox::math::static_vector<Node, Params::kDivsCapacity>;

  // Constructs a window manager on top of ntk using parameters ps and stats st.
  // The stats object is updated on each run.
  WindowManager(Ntk& ntk, Params const& ps, Stats& st)
      : ps_(ps), st_(st), ntk_(ntk), window_(ntk_), times_(ntk) {}

  // Builds a window around node n and updates internal statistics.
  // Returns true iff the resulting window is valid according to ps_.
  bool Run(Node const& n) {
    st_.valid = false;

    if (ntk_.fanin_size(n) == 0) return false;
    if (ntk_.is_pi(n)) return false;

    // Redefine the label associated with the internal part of the window
    window_.Reset(n);

    if (!CollectInps()) return false;
    if (!CollectMffc()) return false;
    if (!CollectDivs()) return false;
    if (!CollectTfos()) return false;
    if (!ExpandInps()) return false;
    if (!FixDominos()) return false;

    SortInps();
    SortDivs();
    SortMffc();
    SortTfos();
    SortOuts();

    FixLabels();

    st_.valid = window_.NumDivs() <= ps_.max_num_divisors;
    st_.valid &= window_.NumInps() <= ps_.kMaxNumLeaves;

    return st_.valid;
  }

  size_t Size() const {
    size_t s = 0;
    s += window_.NumInps();
    s += window_.NumDivs();
    s += window_.NumMffc() * Ntk::max_num_outputs;
    s += window_.NumTfos() * Ntk::max_num_outputs;
    return s;
  }

  // Returns true if the last run produced a valid window.
  bool IsValid() const { return st_.valid; }

  // Returns true if n is an input of the current window.
  bool IsInps(Node const& n) const { return window_.IsInps(n); }

  // Returns true if n is contained anywhere in the current window.
  bool IsContained(Node const& n) const {
    return window_.IsInps(n) || window_.IsDivs(n) || window_.IsMffc(n) ||
           window_.IsTfos(n);
  }
  bool is_contained(Node const& n) const {
    return window_.IsInps(n) || window_.IsDivs(n) || window_.IsMffc(n) ||
           window_.IsTfos(n);
  }

  // Counting helpers.

  size_t NumInps() const { return window_.NumInps(); }
  size_t NumTfos() const { return window_.NumTfos(); }
  size_t NumOuts() const { return window_.NumOuts(); }
  size_t NumMffc() const { return window_.NumMffc(); }
  size_t NumDivs() const { return window_.NumDivs(); }

  void SortInps() { Sort(window_.inps_); }
  void SortDivs() { Sort(window_.divs_); }
  void SortMffc() { Sort(window_.mffc_); }
  void SortTfos() { Sort(window_.tfos_); }
  void SortOuts() { Sort(window_.outs_); }

  template <typename T, uint64_t N>
  void Sort(math::static_vector<T, N>& vec) {
    std::stable_sort(vec.begin(), vec.end(), [&](T const& a, T const& b) {
      return times_.GetTier(a) < times_.GetTier(b);
    });
  }

  // Iterators.

  template <typename Fn>
  void ForeachSignalInps(Fn&& fn) const {
    window_.ForeachSignalInps(fn);
  }

  template <typename Fn>
  void ForeachSignalDivs(Fn&& fn) const {
    window_.ForeachSignalDivs(fn);
  }

  template <typename Fn>
  void ForeachNodeMffc(Fn&& fn) const {
    window_.ForeachNodeMffc(fn);
  }

  template <typename Fn>
  void ForeachNodeTfos(Fn&& fn) const {
    window_.ForeachNodeTfos(fn);
  }

  template <typename Fn>
  void ForeachSignalOuts(Fn&& fn) const {
    window_.ForeachSignalOuts(fn);
  }

  bool CollectInps() {
    FrontNodes front;
    FrontSignals best_inputs;

    size_t num_inps = 0;
    size_t num_divs = 0;

    Node const n = window_.GetRoot();
    window_.MarkSeen(n);

    int32_t tier_last_tfi = times_.GetTier(n) - ps_.max_num_levels_tfi;
    int32_t tier_last_tfo = times_.GetTier(n) + ps_.max_num_levels_tfo;

    ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
      auto const ni = ntk_.get_node(fi);
      if (!window_.WasSeen(ni) && (times_.GetTier(ni) >= tier_last_tfi)) {
        window_.MarkSeen(ni);
        front.push_back(ni);
        num_inps += ntk_.num_outputs(ni);
        num_divs += ntk_.num_outputs(ni);
      }
    });

    best_inputs.clear();
    for (auto const& n : front) {
      ntk_.foreach_output(n, [&](auto const& f) { best_inputs.push_back(f); });
    }

    while (num_inps < (ps_.kMaxNumLeaves + ps_.max_reconvergence)) {
      int best_cost = std::numeric_limits<int>::max();
      std::optional<Node> best_node;

      for (int i = 0; i < static_cast<int>(front.size()); ++i) {
        auto const& n = front[i];
        if (ntk_.is_pi(n)) continue;

        int cost = -static_cast<int>(ntk_.num_outputs(n));
        ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.WasSeen(ni)) cost += ntk_.num_outputs(ni);
          if (times_.GetTier(ni) < tier_last_tfi) cost += ps_.kMaxNumLeaves;
        });

        if (cost < best_cost) {
          best_cost = cost;
          best_node = n;
        }
      }

      if (best_node &&
          ((num_inps + best_cost) <=
           (ps_.kMaxNumLeaves + ps_.max_reconvergence)) &&
          ((num_divs + best_cost + ntk_.num_outputs(*best_node)) <=
           ps_.max_num_divisors)) {
        ntk_.foreach_fanin(*best_node, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.WasSeen(ni)) {
            window_.MarkSeen(ni);
            front.push_back(ni);
            num_inps += ntk_.num_outputs(ni);
            num_divs += ntk_.num_outputs(ni);
          }
        });
        num_inps -= ntk_.num_outputs(*best_node);
        front.erase_value(*best_node);

        if (num_inps <= ps_.kMaxNumLeaves) {
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

    for (auto const& f : best_inputs) window_.MakeInps(ntk_.get_node(f));

    return (window_.NumInps() <= ps_.kMaxNumLeaves) &&
           (window_.divs_.size() <= ps_.max_num_divisors);
  }

  bool CollectMffc() {
    window_.MakeMffc(window_.GetRoot());

    size_t prev_size = 0;
    while (window_.NumMffc() > prev_size) {
      prev_size = window_.NumMffc();

      for (size_t i = 0; i < prev_size; ++i) {
        auto const n = window_.mffc_[i];

        ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);

          // Do not absorb the cut boundary or PIs
          if (window_.IsInps(ni)) return;
          if (ntk_.is_pi(ni)) return;
          if (ntk_.is_po(ni)) return;

          // Standard MFFC test: all fanouts already in the current MFFC
          bool all_in = true;
          ntk_.foreach_fanout(
              ni, [&](auto const& no) { all_in &= window_.IsMffc(no); });

          if (all_in) {
            window_.MakeMffc(ni);
          }
        });
      }

      if (window_.NumMffc() >= Params::kMffcCapacity) break;
    }

    return window_.NumMffc() <= Params::kMffcCapacity;
  }

  bool CollectDivs() {
    if (!CollectTfis()) return false;
    if (!CollectSide()) return false;
    return (window_.NumDivs() <= ps_.max_num_divisors);
  }

  bool CollectTfis() {
    std::stack<Node> stack;

    for (auto const& n : window_.mffc_) {
      ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
        auto const ni = ntk_.get_node(fi);
        if (!window_.WasSeen(ni)) return;
        if (window_.IsInps(ni)) return;
        if (window_.IsMffc(ni)) return;
        stack.push(ni);
      });
    }

    while (!stack.empty()) {
      auto n = stack.top();
      stack.pop();

      if (window_.IsDivs(n) || window_.IsInps(n) || window_.IsMffc(n)) continue;

      window_.MakeDivs(n);

      ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
        auto const ni = ntk_.get_node(fi);
        if (!window_.WasSeen(ni)) return;
        if (window_.IsInps(ni)) return;
        if (window_.IsMffc(ni)) return;
        if (window_.IsDivs(ni)) return;
        stack.push(ni);
      });
    }

    return window_.NumDivs() <= ps_.max_num_divisors;
  }

  bool CollectSide() {
    Node const root = window_.GetRoot();
    int32_t tier_last_tfi = times_.GetTier(root) - ps_.max_num_levels_tfi;
    int32_t tier_last_tfo = times_.GetTier(root) + ps_.max_num_levels_tfo;
    bool expanded = true;
    while (expanded) {
      expanded = false;

      for (int i = 0; i < static_cast<int>(window_.NumDivs()) &&
                      window_.NumDivs() <= ps_.max_num_divisors;
           ++i) {
        auto const& f = window_.divs_[i];
        if (f.output != 0) continue;

        auto const n = ntk_.get_node(f);

        ntk_.foreach_fanout(n, [&](auto const& no) {
          if (times_.GetTier(no) >= times_.GetTier(root)) return;
          if (window_.IsMffc(no)) return;
          if (window_.IsDivs(no)) return;
          if (window_.IsInps(no)) return;

          bool add = true;
          ntk_.foreach_fanin(no, [&](auto const& fi, auto) {
            auto const ni = ntk_.get_node(fi);
            add &= (window_.IsDivs(ni) || window_.IsInps(ni));
          });

          if (add) {
            window_.MakeDivs(no);
            expanded = true;
          }
        });
      }
    }

    return window_.divs_.size() <= ps_.max_num_divisors;
  }

  bool CollectTfos() {
    Node root = window_.GetRoot();
    int32_t tier_last_tfo = times_.GetTier(root) + ps_.max_num_levels_tfo;
    if (ps_.max_num_levels_tfo <= 0) {
      ntk_.foreach_output(
          root, [&](auto const& f) { window_.outs_.insert_unique(f); });
      return true;
    }

    rinox::math::static_vector<Node, Params::kTfosCapacity> current_front;
    rinox::math::static_vector<Node, Params::kTfosCapacity> next_front;

    ntk_.incr_trav_id();
    ntk_.set_visited(root, ntk_.trav_id());

    // -- add immediate fanouts of pivot
    ntk_.foreach_fanout(root, [&](auto const& no) {
      if (ntk_.visited(no) != ntk_.trav_id()) {
        current_front.insert_unique(no);
        ntk_.set_visited(no, ntk_.trav_id());
      }
    });

    // -- count how many NEW inputs are needed
    size_t extra_inputs = 0;
    rinox::math::static_vector<Node, Params::kTfosCapacity> extra_nodes;

    for (auto const& n : current_front) {
      ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
        auto const ni = ntk_.get_node(fi);
        if (!window_.IsInside(ni) && ntk_.visited(ni) != ntk_.trav_id()) {
          extra_inputs += ntk_.num_outputs(ni);
          extra_nodes.push_back(ni);
          ntk_.set_visited(ni, ntk_.trav_id());
        }
      });
    }

    // -- check leaf budget
    if ((window_.NumInps() + extra_inputs) > ps_.kMaxNumLeaves) {
      // Too many inputs: pivot is output instead
      ntk_.foreach_output(
          root, [&](auto const& f) { window_.outs_.insert_unique(f); });
      return true;
    }

    // -- success: commit the nodes
    for (auto const& n : current_front) window_.MakeTfos(n);
    for (auto const& n : extra_nodes) window_.MakeInps(n);

    // -- grow up to num_levels_tfo
    for (int level = 1; level < ps_.max_num_levels_tfo; ++level) {
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
          if (!window_.IsInside(ni) && ntk_.visited(ni) != ntk_.trav_id()) {
            extra_inputs += ntk_.num_outputs(ni);
            extra_nodes.push_back(ni);
            ntk_.set_visited(ni, ntk_.trav_id());
          }
        });
      }

      if ((window_.NumInps() + extra_inputs) > ps_.kMaxNumLeaves) {
        // stop growing here
        break;
      }

      // commit this level
      for (auto const& n : next_front) window_.MakeTfos(n);
      for (auto const& n : extra_nodes) window_.MakeInps(n);

      current_front = next_front;
    }

    if (window_.tfos_.size() > 0) {
      // -- mark outputs: any TFO that leaves window or is a PO
      window_.ForeachNodeTfos([&](auto const& n, auto) {
        bool output = ntk_.is_po(n);
        ntk_.foreach_fanout(n, [&](auto const& no) {
          if (!window_.IsInside(no)) output = true;
        });
        if (output) {
          ntk_.foreach_output(
              n, [&](auto const& f) { window_.outs_.insert_unique(f); });
        }
      });
    } else {
      ntk_.foreach_output(
          root, [&](auto const& f) { window_.outs_.insert_unique(f); });
    }

    return (window_.NumInps() <= ps_.kMaxNumLeaves) &&
           (window_.divs_.size() <= ps_.max_num_divisors) &&
           (window_.tfos_.size() <= ps_.kTfosCapacity) &&
           (window_.outs_.size() <= ps_.kOutsCapacity);
  }

  bool ExpandInps() {
    while (window_.NumInps() < ps_.kMaxNumLeaves) {
      int best_cost = std::numeric_limits<int>::max();
      std::optional<Node> best_node;

      const size_t num_inps = window_.NumInps();
      for (int i = 0; i < static_cast<int>(num_inps); ++i) {
        auto const& f = window_.inps_[i];
        if (f.output != 0) continue;
        auto const n = ntk_.get_node(f);
        if (ntk_.is_pi(n)) continue;

        int cost = -static_cast<int>(ntk_.num_outputs(n));
        ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.IsInps(ni) && !window_.IsDivs(ni))
            cost += ntk_.num_outputs(ni);
        });

        if (cost < best_cost) {
          best_cost = cost;
          best_node = n;
        }
      }

      if (!best_node) break;

      if ((window_.NumInps() + best_cost) <= ps_.kMaxNumLeaves &&
          (window_.NumDivs() + best_cost + ntk_.num_outputs(*best_node)) <=
              ps_.max_num_divisors) {
        window_.MakeDivs(*best_node);
        ntk_.foreach_fanin(*best_node, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.IsInps(ni) && !window_.IsDivs(ni)) window_.MakeInps(ni);
        });
      } else {
        break;
      }
    }

    return (window_.NumInps() <= ps_.kMaxNumLeaves) &&
           (window_.NumDivs() <= ps_.max_num_divisors);
  }

  bool FixDominos() {
    bool any_change = false;
    bool changed = true;

    while (changed) {
      changed = false;

      /* Step 1: refine INPS */
      for (int i = 0; i < static_cast<int>(window_.NumInps()); ++i) {
        auto const f = window_.inps_[i];
        if (f.output != 0) continue;
        auto const n = ntk_.get_node(f);
        if (ntk_.is_pi(n)) continue;

        // Check if all fanins are inside
        bool all_fanins_inside = true;
        ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
          auto const ni = ntk_.get_node(fi);
          if (!window_.IsInps(ni) && !window_.IsDivs(ni))
            all_fanins_inside = false;
        });
        if (!all_fanins_inside) continue;

        // Check if all fanouts are MFFC
        window_.MakeDivs(n);
        changed = true;
        any_change = true;
        break;  // restart since containers have changed
      }
    }

    return (window_.NumInps() <= ps_.kMaxNumLeaves) &&
           (window_.NumDivs() <= ps_.max_num_divisors) &&
           (window_.NumTfos() <= ps_.kTfosCapacity) &&
           (window_.NumOuts() <= ps_.kOutsCapacity);
  }

  void FixLabels() {
    std::set<Signal> inps;
    ForeachSignalInps([&](auto const& f, auto i) {
      window_.MakeInps(ntk_.get_node(f));
      inps.insert(f);
    });
    ForeachSignalDivs([&](auto const& f, auto i) {
      if (inps.find(f) == inps.end()) window_.MakeDivs(ntk_.get_node(f));
    });
  }

  auto const& GetSignalsInps() { return window_.GetSignalsInps(); }
  auto const& GetSignalsDivs() { return window_.GetSignalsDivs(); }
  auto const& GetSignalsOuts() { return window_.GetSignalsOuts(); }
  auto const& GetNodesMffc() { return window_.GetNodesMffc(); }
  auto const& GetNodesTfos() { return window_.GetNodesTfos(); }

  Node GetRoot() { return window_.GetRoot(); }
  Node get_pivot() { return window_.GetRoot(); }
  TimesEvaluator& GetTimes() { return times_; };

 private:
  // Parameters are owned by the caller and must outlive this manager.
  Params const& ps_;
  // Statistics are owned by the caller and updated on each run.
  Stats& st_;
  // Underlying network being windowed.
  Ntk& ntk_;
  // Current window state.
  WindowT window_;
  // Timing evaluator used for node arrival/required times.
  TimesEvaluator times_;
};

}  // namespace windowing
}  // namespace rinox

#endif  // RINOX_WINDOWING_WINDOW_MANAGER_HPP_
