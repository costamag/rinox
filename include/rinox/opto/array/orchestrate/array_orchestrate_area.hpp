/* rinox: C++ logic network library
 * Copyright (C) 20X8-2023  EPFL
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
  \file array_area_resyn.hpp
  \brief Resynthesis using multiple output gates

  \author Andrea Costamagna
*/

#ifndef RINOX_ALGORITHMS_ARRAY_ORCHESTRATE_AREA_HPP_
#define RINOX_ALGORITHMS_ARRAY_ORCHESTRATE_AREA_HPP_

#include <fmt/format.h>

#include <chrono>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>

#include "../../../analyzers/evaluators/levels_evaluator.hpp"
#include "../../../analyzers/trackers/arrival_times_tracker.hpp"
#include "../../../boolean/simd.hpp"
#include "../../../boolean/spfd.hpp"
#include "../../../dependency/window/oracles.hpp"
#include "../../../network/visualize.hpp"
#include "../../../synthesis/lut/synth.hpp"
#include "../../../windowing/simula_manager.hpp"
#include "../../../windowing/window_manager.hpp"
#include "../../../windowing/window_simulator.hpp"
#include "../../profilers/area_profiler.hpp"
#include "../../profilers/delay_profiler.hpp"
#include "../../profilers/power_profiler.hpp"
#include "../../profilers/profilers_utils.hpp"

#ifndef RINOX_NUM_VARS_SIGN
#define RINOX_NUM_VARS_SIGN 6
#endif

#ifndef RINOX_MAX_CUTS_SIZE
#define RINOX_MAX_CUTS_SIZE 6
#endif

#ifndef RINOX_MAX_CUBE_SPFD
#define RINOX_MAX_CUBE_SPFD 12
#endif

#ifndef RINOX_MAX_NUM_LEAVES
#define RINOX_MAX_NUM_LEAVES 12
#endif

namespace rinox {
namespace opto {
namespace array {

struct OrchestrateAreaParamsDefault {
  static constexpr uint32_t kMaxNumLeaves = RINOX_MAX_NUM_LEAVES;
  static constexpr uint32_t kMaxCutsSize = RINOX_MAX_CUTS_SIZE;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
  static constexpr dependency::window::Oracles kSpfdMeth =
      dependency::window::Oracles::kPivs;
  uint32_t stop_after_nth = std::numeric_limits<uint32_t>::max();
  uint32_t max_num_divisors = 1024;
  uint32_t max_reconvergence = 10;
  int32_t max_num_levels_tfo = 0;
  bool try_0Bresub = false;
  bool try_0Xresub = false;
  bool try_1Bresub = false;
  bool try_1Xresub = false;
  bool preserve_time = false;
};

struct OrchestrateAreaStatsDefault {
  mockturtle::stopwatch<>::duration time_window_construction{0};
  mockturtle::stopwatch<>::duration time_window_simulation{0};
  mockturtle::stopwatch<>::duration time_0Bresub{0};
  mockturtle::stopwatch<>::duration time_0Xresub{0};
  mockturtle::stopwatch<>::duration time_1Bresub{0};
  mockturtle::stopwatch<>::duration time_1Xresub{0};
  uint32_t num_changes{0};
  uint32_t num_0Bresub{0};
  uint32_t num_0Xresub{0};
  uint32_t num_1Bresub{0};
  uint32_t num_1Xresub{0};

  double GetTimeWindowConstr() const {
    return mockturtle::to_seconds(time_window_construction);
  }
  double GetTimeWindowSimula() const {
    return mockturtle::to_seconds(time_window_simulation);
  }
  double GetTime0BResub() const { return mockturtle::to_seconds(time_0Bresub); }
  double GetTime0XResub() const { return mockturtle::to_seconds(time_0Xresub); }
  double GetTime1BResub() const { return mockturtle::to_seconds(time_1Bresub); }
  double GetTime1XResub() const { return mockturtle::to_seconds(time_1Xresub); }

  double GetTimeTotal() const {
    return GetTimeWindowConstr() + GetTimeWindowSimula() + GetTimeOpto();
  }

  double GetTimeOpto() const {
    return GetTime0BResub() + GetTime0XResub() + GetTime1BResub() +
           GetTime1XResub();
  }

  void Record0BResub() { num_0Bresub++; }
  void Record0XResub() { num_0Xresub++; }
  void Record1BResub() { num_1Bresub++; }
  void Record1XResub() { num_1Xresub++; }

  uint32_t GetNum0BResub() const { return num_0Bresub; }
  uint32_t GetNum0XResub() const { return num_0Xresub; }
  uint32_t GetNum1BResub() const { return num_1Bresub; }
  uint32_t GetNum1XResub() const { return num_1Xresub; }

  uint32_t GetNumOptos() const {
    return GetNum0BResub() + GetNum0XResub() + GetNum1BResub() +
           GetNum1XResub();
  }

  void Report() const {
    std::cout << fmt::format("================ TIME  ====================\n");
    std::cout << fmt::format("    time window construction       = {:5.2f}\n",
                             GetTimeWindowConstr());
    std::cout << fmt::format("    time window simulation         = {:5.2f}\n",
                             GetTimeWindowSimula());

    std::cout << fmt::format("    time 0-cost basic version      = {:5.2f}\n",
                             GetTime0BResub());
    std::cout << fmt::format("    time 0-cost with augmentation  = {:5.2f}\n",
                             GetTime0XResub());
    std::cout << fmt::format("    time 1-cost basic version      = {:5.2f}\n",
                             GetTime1BResub());
    std::cout << fmt::format("    time 1-cost with augmentation  = {:5.2f}\n",
                             GetTime1XResub());

    std::cout << fmt::format("    time total                     = {:5.2f}\n",
                             GetTimeTotal());
    std::cout << fmt::format("================ QUALITY ==================\n");
    std::cout << fmt::format("    count 0-cost basic version     = {:5d}\n",
                             GetNum0BResub());
    std::cout << fmt::format("    count 0-cost with augmentation = {:5d}\n",
                             GetNum0XResub());
    std::cout << fmt::format("    count 1-cost basic version     = {:5d}\n",
                             GetNum1BResub());
    std::cout << fmt::format("    count 1-cost with augmentation = {:5d}\n",
                             GetNum1XResub());
    std::cout << fmt::format("    count total                    = {:5d}\n",
                             GetNumOptos());
  }
};

namespace detail {

template <class Ntk, typename Params, typename Stats>
class OrchestrateAreaImpl {
  struct WindowManagerParams {
    uint32_t max_num_divisors = 1024;
    uint32_t max_reconvergence = 10;
    int32_t max_num_levels_tfo = 0;
    static constexpr uint32_t max_num_leaves = Params::kMaxNumLeaves;
  };

  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using Signature = kitty::static_truth_table<Params::kMaxNumLeaves>;
  using WindowManager = windowing::window_manager<Ntk, WindowManagerParams>;
  using WindowSimulator =
      windowing::window_simulator<Ntk, Params::kMaxNumLeaves>;
  using Spfd = boolean::spfd<Signature, 1u << Params::kMaxCutsSize>;
  static constexpr uint32_t kMaxNumOutputs = Ntk::kMaxNumOutputs;
  using Dependencies = typename dependency::window::OracleDispatch<
      Params::kSpfdMeth, Ntk, WindowManager, WindowSimulator,
      Params::kMaxCutsSize>::type;

 public:
  OrchestrateAreaImpl(Ntk& ntk, Params ps, Stats& st)
      :  // parameters
        ps_(ps),
        window_manager_ps_{ps_.max_num_divisors, ps_.max_reconvergence,
                           ps_.max_num_levels_tfo},
        // statistics
        st_(st),
        // attributes
        ntk_(ntk),
        window_manager_(ntk_, window_manager_ps_, window_manager_st_),
        window_simulator_(ntk_),
        dependencies_(ntk, window_manager_, window_simulator_),
        levels_(ntk) {}

  void Report();
  void Run();

 private:
  bool SkipNode(Node) const;
  bool Try0BResub(Node);
  bool Try0XResub(Node);
  bool Try1BResub(Node);
  bool Try1XResub(Node);
  bool Try1XResubSolo(Node);
  bool Try1XResubDual(Node);

  std::optional<Node> FindBestAugmentableNode(Spfd& spfd);
  std::optional<Signal> AugmentNode(Node, Signature const&);
  std::optional<Signal> AugmentNode(Node, Signature const&, Signature const&);
  kitty::dynamic_truth_table Interpolate(std::vector<Signal> const&,
                                         Signature const&);
  kitty::dynamic_truth_table Interpolate(std::vector<Signal> const&,
                                         Signature const&, Signature const&);

 private:
  // Parameters
  Params ps_;
  WindowManagerParams window_manager_ps_;
  // Statistics
  Stats& st_;
  windowing::window_manager_stats window_manager_st_;
  windowing::SimulaManagerStats window_simulator_st_;
  // Attributes
  Ntk& ntk_;
  Spfd spfd_;
  WindowManager window_manager_;
  WindowSimulator window_simulator_;
  Dependencies dependencies_;
  rinox::analyzers::evaluators::LevelsEvaluator<Ntk> levels_;
  rinox::boolean::Smoother<Signature, int, Params::kMaxInpsSolo> smoother_;
  rinox::math::SpinSystem<Signature, 1u << Params::kMaxInpsSolo, kMaxNumOutputs>
      spins_;
};

template <typename Ntk, typename Params, typename Stats>
void OrchestrateAreaImpl<Ntk, Params, Stats>::Report() {
  std::vector<uint32_t> mffc_sizes_solo;
  std::vector<uint32_t> mffc_sizes_dual;
  ntk_.foreach_gate([&](auto const& n) {
    window_manager_.run(n);
    uint32_t mffc_size = window_manager_.num_mffc();
    if (ntk_.num_outputs(n) == 1) {
      while (mffc_sizes_solo.size() <= mffc_size) mffc_sizes_solo.push_back(0);
      mffc_sizes_solo[mffc_size] += 1;
    } else {
      while (mffc_sizes_dual.size() <= mffc_size) mffc_sizes_dual.push_back(0);
      mffc_sizes_dual[mffc_size] += 1;
    }
    return true;
  });
  for (uint32_t s = 0; s < mffc_sizes_solo.size(); ++s) {
    std::cout << s << " 1 : " << mffc_sizes_solo[s] << std::endl;
  }
  for (uint32_t s = 0; s < mffc_sizes_dual.size(); ++s) {
    std::cout << s << " 2 : " << mffc_sizes_dual[s] << std::endl;
  }
}

// Implementation of Run split into steps
template <typename Ntk, typename Params, typename Stats>
void OrchestrateAreaImpl<Ntk, Params, Stats>::Run() {
  levels_.ComputeLevels();
  int const last_level = levels_.GetLast();

  for (int level = last_level; level > 0; --level) {
    bool update = false;
    levels_.ForeachNode(level, [&](auto const& n) {
      if (SkipNode(n)) return true;

      mockturtle::call_with_stopwatch(st_.time_window_construction,
                                      [&]() { window_manager_.run(n); });
      mockturtle::call_with_stopwatch(st_.time_window_simulation, [&]() {
        window_simulator_.run(window_manager_);
      });

      if (!window_manager_.is_valid()) return true;

      if (ps_.try_0Bresub &&
          mockturtle::call_with_stopwatch(st_.time_0Bresub,
                                          [&]() { return Try0BResub(n); })) {
        st_.Record0BResub();
        update = true;
        return true;
      }
      if (ps_.try_0Xresub &&
          mockturtle::call_with_stopwatch(st_.time_0Xresub,
                                          [&]() { return Try0XResub(n); })) {
        st_.Record0XResub();
        update = true;
        return true;
      }
      if (ps_.try_1Bresub &&
          mockturtle::call_with_stopwatch(st_.time_1Bresub,
                                          [&]() { return Try1BResub(n); })) {
        st_.Record1BResub();
        update = true;
        return true;
      }
      if (ps_.try_1Xresub &&
          mockturtle::call_with_stopwatch(st_.time_1Xresub,
                                          [&]() { return Try1XResub(n); })) {
        st_.Record1XResub();
        update = true;
        return true;
      }
      return true;
    });

    // update the levels if there was an optimization
    if (update) levels_.ComputeLevels();
  }
}

template <typename Ntk, typename Params, typename Stats>
bool OrchestrateAreaImpl<Ntk, Params, Stats>::SkipNode(Node n) const {
  if (ntk_.fanout_size(n) <= 0) return true;
  if (ntk_.is_pi(n) || ntk_.is_constant(n)) return true;
  if (ntk_.is_dead(n)) return true;
  return false;
}

template <typename Ntk, typename Params, typename Stats>
bool OrchestrateAreaImpl<Ntk, Params, Stats>::Try0BResub(Node n) {
  std::vector<Signal> new_signals;
  ntk_.foreach_output(n, [&](auto const& f) {
    if (new_signals.size() != f.output) return;

    int tA = std::numeric_limits<int>::max();
    auto const& onset = window_simulator_.Get(f);
    std::optional<Signal> resub;
    window_manager_.foreach_divisor([&](auto const& d, auto) {
      if (!resub || (levels_.GetArrival(d) < tA)) {
        auto const& divisor = window_simulator_.Get(d);
        if (kitty::equal(onset, divisor)) {
          resub = d;
        }
      }
    });
    if (resub) new_signals.push_back(*resub);
  });

  if (new_signals.size() == ntk_.num_outputs(n)) {
    ntk_.substitute_node(n, new_signals);
    return true;
  }
  return false;
}

template <typename Ntk, typename Params, typename Stats>
bool OrchestrateAreaImpl<Ntk, Params, Stats>::Try0XResub(Node n) {
  if (ntk_.num_outputs(n) > 1) return false;
  std::vector<Signal> new_children;
  kitty::dynamic_truth_table new_function;

  Signature target = window_simulator_.Get(ntk_.make_signal(n));

  spfd_.init(target);
  auto node = FindBestAugmentableNode(spfd_);
  if (node) {
    auto const fnew = AugmentNode(*node, target);
    if (fnew) {
      Node nnew = ntk_.get_node(*fnew);
      window_simulator_.compute(window_manager_, nnew);
      bool cec1 = kitty::equal(window_simulator_.Get(Signal{nnew, 0}),
                               window_simulator_.Get(ntk_.make_signal(*node)));
      bool cec2 = kitty::equal(window_simulator_.Get(Signal{nnew, 1}), target);

      if (cec1 && cec2) {
        ntk_.substitute_node(*node, Signal{nnew, 0});
        ntk_.substitute_node(n, Signal{nnew, 1});
        return true;
      } else {
        std::cout << "Failed cec 0XResub " << std::endl;
        ntk_.take_out_node(nnew);
      }
    }
  }
  return false;
}

template <typename Ntk, typename Params, typename Stats>
bool OrchestrateAreaImpl<Ntk, Params, Stats>::Try1BResub(Node n) {
  if (window_manager_.num_mffc() <= 1) return false;

  dependencies_.Init();
  size_t cut_size =
      ntk_.num_outputs(n) == 1 ? Params::kMaxInpsSolo : Params::kMaxInpsDual;
  auto cut = dependencies_.GetNext(cut_size);
  if (cut) {
    std::vector<Signal> children;
    for (auto f : *cut) {
      children.push_back(f);
    }
    if (ntk_.num_outputs(n) == 1) {
      auto const& target = window_simulator_.Get(ntk_.make_signal(n));
      auto tt = Interpolate(children, target);
      auto fnew = ntk_.create_node(children, tt);
      auto nnew = ntk_.get_node(fnew);
      window_simulator_.compute(window_manager_, nnew);
      if (kitty::equal(window_simulator_.Get(fnew), target)) {
        ntk_.substitute_node(n, fnew);
        return true;
      } else {
        ntk_.take_out_node(nnew);
        std::cout << "Failed cec 1 Resub " << std::endl;
      }
    } else if (ntk_.num_outputs(n) == 2) {
      auto const& target0 = window_simulator_.Get(Signal{n, 0});
      auto const& target1 = window_simulator_.Get(Signal{n, 1});
      auto tt0 = Interpolate(children, target0);
      auto tt1 = Interpolate(children, target1);
      auto fnew = ntk_.create_node(children, {tt0, tt1});
      auto nnew = ntk_.get_node(fnew);
      window_simulator_.compute(window_manager_, nnew);
      bool cec1 = kitty::equal(target0, window_simulator_.Get(Signal{nnew, 0}));
      bool cec2 = kitty::equal(target1, window_simulator_.Get(Signal{nnew, 1}));

      if (cec1 && cec2) {
        std::vector<Signal> new_signals = {Signal{nnew, 0}, Signal{nnew, 1}};
        ntk_.substitute_node(n, new_signals);
        return true;
      } else {
        std::cout << "Failed cec 1 Resub dual " << std::endl;
        ntk_.take_out_node(nnew);
      }
    }
  }
  return false;
}

template <typename Ntk, typename Params, typename Stats>
bool OrchestrateAreaImpl<Ntk, Params, Stats>::Try1XResub(Node n) {
  if (window_manager_.num_mffc() <= 1) return false;
  if (ntk_.num_outputs(n) == 1)
    return Try1XResubSolo(n);
  else if (ntk_.num_outputs(n) == 2)
    return Try1XResubDual(n);
  return false;
}

template <typename Ntk, typename Params, typename Stats>
bool OrchestrateAreaImpl<Ntk, Params, Stats>::Try1XResubSolo(Node n) {
  if (window_manager_.num_mffc() <= 1) return false;
  Signature target = window_simulator_.Get(ntk_.make_signal(n));
  Signature func = window_simulator_.Get(ntk_.make_signal(n));
  Signature care = ~func.construct();

  std::vector<Node> added;
  std::vector<std::pair<Node, Signal>> todos;
  std::vector<Signal> new_fanin;

  ntk_.incr_trav_id();
  spfd_.init(func);
  std::set<Node> avoid;
  int cnt = 0;
  while ((cnt++ < Params::kMaxInpsSolo) &&
         (new_fanin.size() < Params::kMaxInpsSolo) && !(spfd_.is_covered())) {
    std::optional<Signal> best_asis;
    std::optional<Node> best_dupl;
    uint32_t best_cost_asis = std::numeric_limits<uint32_t>::max();
    uint32_t best_cost_dupl = std::numeric_limits<uint32_t>::max();
    window_manager_.foreach_divisor([&](auto const& fd, auto) {
      auto const nd = ntk_.get_node(fd);
      if ((avoid.find(nd) == avoid.end()) &&
          (levels_.GetArrival(fd) < levels_.GetArrival(ntk_.make_signal(n)))) {
        uint32_t const cost = spfd_.evaluate(window_simulator_.Get(fd));
        if (cost < best_cost_asis) {
          best_cost_asis = cost;
          best_asis = fd;
        }
        if ((!window_manager_.is_input(nd)) && (ntk_.num_outputs(nd) == 1) &&
            (cost < best_cost_dupl) &&
            (ntk_.fanin_size(nd) <= Params::kMaxInpsDual)) {
          best_cost_dupl = cost;
          best_dupl = nd;
        }
      }
    });
    if (best_dupl) {
      auto fnew = AugmentNode(*best_dupl, func, care);
      avoid.insert(*best_dupl);
      if (fnew) {
        Node nnew = ntk_.get_node(*fnew);
        window_simulator_.compute(window_manager_, nnew);
        auto cost = spfd_.evaluate(window_simulator_.Get(Signal{nnew, 1}));

        if (cost < best_cost_asis) {
          added.push_back(nnew);
          todos.push_back(std::pair{*best_dupl, Signal{nnew, 0}});
          new_fanin.push_back(Signal{nnew, 1});
          ntk_.set_visited(*best_dupl, ntk_.trav_id());
        } else if (best_asis) {
          if (ntk_.fanout_size(nnew) == 0) ntk_.take_out_node(nnew);

          new_fanin.push_back(*best_asis);
        } else {
          return false;
        }
      } else if (best_asis) {
        new_fanin.push_back(*best_asis);
      }
    } else if (best_asis) {
      new_fanin.push_back(*best_asis);
    } else {
      return false;
    }

    spfd_.init(target);
    for (auto f : new_fanin) {
      spfd_.update(window_simulator_.Get(f));
    }
    if (spfd_.is_covered()) break;

    spins_.clear();
    care ^= care;
    ;
    for (uint32_t m = 0; m < spfd_.num_masks_; ++m) {
      if (!spfd_.is_killed(m)) {
        spins_.AddSpin(target, spfd_.masks_[m], m);
        care |= spfd_.masks_[m];
      }
    }
    spins_.MinimizeEnergy();
    func = target;
    spins_.ForeachSpin([&](auto const& s) {
      if (s.value < 0) func ^= s.care;
    });
    func = smoother_.RunLabelPropagation(func, care);
  }

  if (spfd_.is_covered()) {
    auto tt =
        Interpolate(new_fanin, window_simulator_.Get(ntk_.make_signal(n)));
    auto fnew = ntk_.create_node(new_fanin, tt);
    window_simulator_.compute(window_manager_, fnew.index);
    if (!kitty::equal(window_simulator_.Get(ntk_.make_signal(n)),
                      window_simulator_.Get(fnew))) {
      std::cout << "Failed cec 1 Rexub " << std::endl;
      ntk_.take_out_node(fnew.index);
      return false;
    }
    ntk_.substitute_node(n, fnew);
    for (auto const p : todos) ntk_.substitute_node(p.first, p.second);
    return true;
  }
  for (auto nd : added) {
    if (ntk_.fanout_size(nd) == 0) ntk_.take_out_node(nd);
  }
  return false;
}

template <typename Ntk, typename Params, typename Stats>
bool OrchestrateAreaImpl<Ntk, Params, Stats>::Try1XResubDual(Node n) {
  return false;
}

template <typename Ntk, typename Params, typename Stats>
auto OrchestrateAreaImpl<Ntk, Params, Stats>::FindBestAugmentableNode(
    Spfd& spfd) -> std::optional<Node> {
  std::optional<Node> res;
  uint32_t best_info_cost = std::numeric_limits<uint32_t>::max();
  uint32_t best_size_cost = Params::kMaxInpsDual;
  int best_time_cost = std::numeric_limits<int>::max();

  window_manager_.foreach_divisor([&](auto const& fd, auto) {
    auto nd = ntk_.get_node(fd);
    uint32_t const size_cost = ntk_.fanin_size(nd);
    int const time_cost = levels_.GetArrival(fd);

    if (!window_manager_.is_input(nd) && (ntk_.num_outputs(nd) == 1) &&
        (size_cost <= best_size_cost) && (time_cost <= best_time_cost)) {
      uint32_t const info_cost = spfd.evaluate(window_simulator_.Get(fd));
      if (info_cost < best_info_cost) {
        best_info_cost = info_cost;
        best_time_cost = time_cost;
        best_size_cost = size_cost;
        res = nd;
      }
    }
  });
  return res;
}

template <typename Ntk, typename Params, typename Stats>
auto OrchestrateAreaImpl<Ntk, Params, Stats>::AugmentNode(
    Node n, Signature const& target) -> std::optional<Signal> {
  Signature care = target | ~target;
  return AugmentNode(n, target, care);
}

template <typename Ntk, typename Params, typename Stats>
auto OrchestrateAreaImpl<Ntk, Params, Stats>::AugmentNode(
    Node n, Signature const& target,
    Signature const& care) -> std::optional<Signal> {
  auto f = ntk_.make_signal(n);

  auto children = ntk_.get_children(n);
  spfd_.init(target, care);
  for (auto fi : children) spfd_.update(window_simulator_.Get(fi));

  if (children.size() < Params::kMaxInpsDual) {
    ntk_.incr_trav_id();
    window_manager_.MarkTfoInside(n, ntk_.trav_id());

    while ((children.size() < Params::kMaxInpsDual) && !spfd_.is_covered()) {
      std::optional<Signal> best_div;
      uint32_t best_cost = std::numeric_limits<uint32_t>::max();
      window_manager_.foreach_divisor([&](auto const& fd, auto) {
        auto nd = ntk_.get_node(fd);
        if ((ntk_.visited(nd) != ntk_.trav_id()) &&
            ((levels_.GetArrival(fd) <
              levels_.GetArrival(f))))  //|| !ExistsPath(n,nd)))
        {
          uint32_t const cost = spfd_.evaluate(window_simulator_.Get(fd));
          if (cost < best_cost) {
            best_cost = cost;
            best_div = fd;
          }
        }
      });
      if (best_div) {
        children.push_back(*best_div);
        spfd_.update(window_simulator_.Get(*best_div));
      } else
        return std::nullopt;
    }
  }

  if (spfd_.is_covered()) {
    auto tt1 =
        Interpolate(children, window_simulator_.Get(ntk_.make_signal(n)));
    auto tt2 = Interpolate(children, target, care);
    auto fnew = ntk_.create_node(children, {tt1, tt2});
    return fnew;
  }
  return std::nullopt;
}

template <typename Ntk, typename Params, typename Stats>
auto OrchestrateAreaImpl<Ntk, Params, Stats>::Interpolate(
    std::vector<Signal> const& children,
    Signature const& target) -> kitty::dynamic_truth_table {
  Signature care = target | ~target;
  return Interpolate(children, target, care);
}

template <typename Ntk, typename Params, typename Stats>
auto OrchestrateAreaImpl<Ntk, Params, Stats>::Interpolate(
    std::vector<Signal> const& children, Signature const& target,
    Signature const& care) -> kitty::dynamic_truth_table {
  uint32_t const num_vars = children.size();
  kitty::dynamic_truth_table tt(num_vars);
  uint32_t const num_minterms = 1u << num_vars;
  Signature mask;
  for (uint32_t m = 0; m < num_minterms; ++m) {
    mask |= ~mask;
    for (uint32_t v = 0; v < num_vars; ++v) {
      if (((m >> v) & 0x1) > 0)
        mask &= window_simulator_.Get(children[v]);
      else
        mask &= ~window_simulator_.Get(children[v]);
    }
    uint32_t c1 = kitty::count_ones(care & mask & target);
    uint32_t c0 = kitty::count_ones(care & mask & ~target);
    if (c0 == 0 && c1 > 0)
      kitty::set_bit(tt, m);
    else if (c1 == 0)
      kitty::clear_bit(tt, m);
    else {
      std::cout << "Warning: don't know how to interpolate" << std::endl;
    }
  }
  return tt;
}

}  // namespace detail

// Public API
template <class Ntk, typename Params = OrchestrateAreaParamsDefault,
          typename Stats = OrchestrateAreaStatsDefault>
void OrchestrateArea(Ntk& ntk, Params ps = {}, Stats* pst = nullptr) {
  Stats st;
  detail::OrchestrateAreaImpl<Ntk, Params, Stats> p(ntk, ps, st);
  p.Run();
  st.Report();
  if (pst) *pst = st;
}

}  // namespace array
}  // namespace opto
}  // namespace rinox

#endif  // RINOX_ALGORITHMS_ARRAY_AREA_ORCHESTRA_HPP_
