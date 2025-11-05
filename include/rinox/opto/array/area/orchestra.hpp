/* rinox: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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

#ifndef RINOX_ALGORITHMS_ARRAY_AREA_ORCHESTRA_HPP_
#define RINOX_ALGORITHMS_ARRAY_AREA_ORCHESTRA_HPP_

#include <fmt/format.h>

#include <chrono>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <optional>

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

struct DefaultAreaResynStaticParams {
  static constexpr uint32_t kMaxNumLeavesXL = 20u;
  static constexpr uint32_t kMaxNumLeaves = RINOX_MAX_NUM_LEAVES;
  static constexpr uint32_t kNumVarsSign = RINOX_NUM_VARS_SIGN;
  static constexpr uint32_t kMaxCutsSize = RINOX_MAX_CUTS_SIZE;
  static constexpr uint32_t kMaxCubeSpfd = RINOX_MAX_CUBE_SPFD;
  static constexpr uint32_t kMaxInpsDual = 5u;
  static constexpr uint32_t kMaxInpsSolo = 6u;
  static constexpr synthesis::lut::LutSynthMethod kLutSynth =
      synthesis::lut::LutSynthMethod::kShannonSolo;
  static constexpr dependency::window::Oracles kSpfdMeth =
      dependency::window::Oracles::kPivs;
};

struct AreaResynParams {
  uint32_t stop_after_nth = std::numeric_limits<uint32_t>::max();
  uint32_t max_num_divisors = 128;
  uint32_t max_reconvergence = 10;
  int32_t odc_levels = 0;
  bool cec = false;
  bool try_0resub_solo = false;
  bool try_1resub_solo = false;
  bool try_2resub_solo = false;
  bool try_1resub_dual = false;
  bool try_1rexub_solo = false;
  bool try_0rexub_solo = false;
  bool try_merger_solo = false;
};

struct AreaResynStats {
  mockturtle::stopwatch<>::duration time_window_constr{0};
  mockturtle::stopwatch<>::duration time_window_simula{0};
  mockturtle::stopwatch<>::duration time_0resub_solo{0};
  mockturtle::stopwatch<>::duration time_1resub_solo{0};
  mockturtle::stopwatch<>::duration time_2resub_solo{0};
  mockturtle::stopwatch<>::duration time_1resub_dual{0};
  mockturtle::stopwatch<>::duration time_1rexub_solo{0};
  mockturtle::stopwatch<>::duration time_0rexub_solo{0};
  mockturtle::stopwatch<>::duration time_merger_solo{0};
  uint32_t num_changes{0};
  uint32_t num_0resub_solo{0};
  uint32_t num_1resub_solo{0};
  uint32_t num_2resub_solo{0};
  uint32_t num_1resub_dual{0};
  uint32_t num_1rexub_solo{0};
  uint32_t num_0rexub_solo{0};
  uint32_t num_merger_solo{0};
  synthesis::lut::LutSynthStats st_synth;

  double GetTimeWindowConstr() const {
    return mockturtle::to_seconds(time_window_constr);
  }
  double GetTimeWindowSimula() const {
    return mockturtle::to_seconds(time_window_simula);
  }
  double GetTimeMergerSolo() const {return mockturtle::to_seconds(time_merger_solo);}

  double GetTime0ResubSolo() const {return mockturtle::to_seconds(time_0resub_solo);}
  double GetTime1ResubSolo() const {return mockturtle::to_seconds(time_1resub_solo);}
  double GetTime2ResubSolo() const {return mockturtle::to_seconds(time_2resub_solo);}

  double GetTime0RexubSolo() const {return mockturtle::to_seconds(time_0rexub_solo);}
  double GetTime1RexubSolo() const {return mockturtle::to_seconds(time_1rexub_solo);}

  double GetTime1ResubDual() const {return mockturtle::to_seconds(time_1resub_dual);}

  double GetTimeTotal() const {
    return GetTimeWindowConstr() + GetTimeWindowSimula() + GetTimeOptoSolo();
  }

  double GetTimeOptoSolo() const {
    return GetTime0ResubSolo() + GetTime1ResubSolo() + GetTime2ResubSolo() + 
           GetTime1ResubDual() + 
           GetTime0RexubSolo() + GetTime1RexubSolo() +
           GetTimeMergerSolo();
  }

  void Record0ResubSolo() { num_0resub_solo++; }
  void Record1ResubSolo() { num_1resub_solo++; }
  void Record2ResubSolo() { num_2resub_solo++; }

  void Record0RexubSolo() { num_0rexub_solo++; }
  void Record1RexubSolo() { num_1rexub_solo++; }

  void Record1ResubDual() { num_1resub_dual++; }

  void RecordMergerSolo() { num_merger_solo++; }

  uint32_t GetNum0ResubSolo() const { return num_0resub_solo; }
  uint32_t GetNum1ResubSolo() const { return num_1resub_solo; }
  uint32_t GetNum2ResubSolo() const { return num_2resub_solo; }

  uint32_t GetNum0RexubSolo() const { return num_0rexub_solo; }
  uint32_t GetNum1RexubSolo() const { return num_1rexub_solo; }

  uint32_t GetNum1ResubDual() const { return num_1resub_dual; }

  uint32_t GetNumMergerSolo() const { return num_merger_solo; }

  uint32_t GetNumOptos() const {
    return GetNum0ResubSolo() + GetNum1ResubSolo() + GetNum2ResubSolo() + 
           GetNum0RexubSolo() + GetNum1RexubSolo() +
           GetNum1ResubDual() +
           GetNumMergerSolo();
  }

  void Report() const {
    std::cout << fmt::format("================ TIME  ====================\n");
    std::cout << fmt::format("    time window construction      = {:5.2f}\n",
                             GetTimeWindowConstr());
    std::cout << fmt::format("    time window simulation        = {:5.2f}\n",
                             GetTimeWindowSimula());

    std::cout << fmt::format("    time 0-resub solo             = {:5.2f}\n",
                             GetTime0ResubSolo());
    std::cout << fmt::format("    time 1-resub solo             = {:5.2f}\n",
                             GetTime1ResubSolo());
    std::cout << fmt::format("    time 2-resub solo             = {:5.2f}\n",
                             GetTime2ResubSolo());

    std::cout << fmt::format("    time 0-rexub solo             = {:5.2f}\n",
                             GetTime0RexubSolo());
    std::cout << fmt::format("    time 1-rexub solo             = {:5.2f}\n",
                            GetTime1RexubSolo());

    std::cout << fmt::format("    time merger solo              = {:5.2f}\n",
                             GetTimeMergerSolo());

    std::cout << fmt::format("    time 1-resub dual             = {:5.2f}\n",
                             GetTime1ResubDual());

    std::cout << fmt::format("    time total                    = {:5.2f}\n",
                             GetTimeTotal());
    std::cout << fmt::format("================ QUALITY ==================\n");
    std::cout << fmt::format("    count 0-resub solo            = {:5d}\n",
                             GetNum0ResubSolo());
    std::cout << fmt::format("    count 1-resub solo            = {:5d}\n",
                             GetNum1ResubSolo());
    std::cout << fmt::format("    count 2-resub solo            = {:5d}\n",
                             GetNum2ResubSolo());

    std::cout << fmt::format("    count 0-rexub solo            = {:5d}\n",
                             GetNum0RexubSolo());
    std::cout << fmt::format("    count 1-rexub solo            = {:5d}\n",
                             GetNum1RexubSolo());

    std::cout << fmt::format("    count 1-resub dual            = {:5d}\n",
                             GetNum1ResubDual());

    std::cout << fmt::format("    count merger solo             = {:5d}\n",
                             GetNumMergerSolo());

    std::cout << fmt::format("    count total                   = {:5d}\n",
                             GetNumOptos());
    st_synth.Report();
  }
};

template <typename Ntk>
struct ReturnOptoData {
  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;

  Node old_node;
  Signal new_signal;
};

namespace detail {

template <class Ntk, typename StaticParams = DefaultAreaResynStaticParams>
class AreaResynImpl {
  struct WindowManagerParams {
    static constexpr uint32_t max_num_leaves = StaticParams::kMaxNumLeaves;
    uint32_t max_num_divisors;
    uint32_t max_reconvergence;
    int32_t odc_levels;
  };

  struct SimulaManagerParams {
    uint32_t max_num_divisors = 128;
    uint32_t max_reconvergence = 10;
    int32_t odc_levels = 0;
    static constexpr uint32_t kMaxNumLeaves = StaticParams::kMaxNumLeavesXL;
    static constexpr uint32_t kNumVarsSigns = StaticParams::kNumVarsSign;
  };

  using Node = typename Ntk::node;
  using Signal = typename Ntk::signal;
  using ArrivalsTracker = analyzers::trackers::arrival_times_tracker<Ntk>;
  using Signature = kitty::static_truth_table<StaticParams::kMaxNumLeaves>;
  using SimulaManager = windowing::SimulaManager<Ntk, SimulaManagerParams>;
  using WindowManager = windowing::window_manager<Ntk, WindowManagerParams>;
  using WindowSimulator =
      windowing::window_simulator<Ntk, StaticParams::kMaxNumLeaves>;
  using Spfds =
      boolean::spfd<Signature, 1u << StaticParams::kMaxCutsSize>;

  static constexpr uint32_t kMaxNumOutputs = Ntk::kMaxNumOutputs;

  struct SharedOptoWindow {
    WindowManager manager;
    WindowSimulator simulator;
    SimulaManager simula;
    Node root{};
  };

  struct LutSynthParams : synthesis::lut::DefaultLutSynthParams {
    static constexpr uint32_t kNumVarsCube = StaticParams::kMaxNumLeaves;
    static constexpr uint32_t kMaxInpsSolo = StaticParams::kMaxInpsSolo;
    static constexpr uint32_t kMaxInpsDual = StaticParams::kMaxInpsDual;
    static constexpr uint32_t kMaxNumGates = 32u;
  };

  using LutSynth =
      typename synthesis::lut::LutSynthDispatch<StaticParams::kLutSynth, Ntk,
                                                WindowManager, WindowSimulator,
                                                LutSynthParams>::type;

  using SpfdMeth = typename dependency::window::OracleDispatch<
      StaticParams::kSpfdMeth, Ntk, WindowManager, WindowSimulator,
      StaticParams::kMaxCutsSize>::type;

 public:
  AreaResynImpl(Ntk& ntk, AreaResynParams ps, AreaResynStats& st)
      : ntk_(ntk),
        ps_(ps),
        st_(st),
        wm_params_{ps_.max_num_divisors, ps_.max_reconvergence, ps_.odc_levels},
        sm_params_{ps_.max_num_divisors, ps_.max_reconvergence, ps_.odc_levels},
        window_{WindowManager(ntk_, wm_params_, wm_stats_),
                WindowSimulator(ntk_),
                SimulaManager(ntk_, sm_params_, sm_stats_)},
        synth_(ntk, window_.manager, window_.simulator, st_.st_synth),
        spfds_(ntk, window_.manager, window_.simulator),
        arrivals_(ntk_) {
    for (int n = 0; n <= StaticParams::kMaxInpsSolo; ++n) {
      proj_fns_.emplace_back();
      for (int i = 0; i < n; ++i) {
        proj_fns_[n].emplace_back(n);
        kitty::create_nth_var(proj_fns_[n][i], i);
      }
    }
  }

  void Run();
  void ReportNtkState();

 private:
  bool SkipNode(Node) const;
  bool ExistsPath(Node, Node) const;
  bool NeedWindow() const;

  void OptimizeSolo(Node);

  bool Try0ResubSolo(Node);
  bool Try1ResubSolo(Node);
  bool Try2ResubSolo(Node);

  bool Try0RexubSolo(Node);
  bool Try1RexubSolo(Node);

  bool Try1ResubDual(Node);
  bool TryMergerSolo(Node);

  std::optional<Signal> AddOutput(Node, Signature const&);
  std::optional<Signal> AddOutput(Node, Signature const&, Signature const&);
  kitty::dynamic_truth_table Interpolate(std::vector<Signal> const&,
                                         Signature const&, Signature const&);
  kitty::dynamic_truth_table Interpolate(std::vector<Signal> const&,
                                         Signature const&);

  void OptimizeDual(Node);
  void Merge(Node, Node);
  kitty::dynamic_truth_table Approximate(Signature const&, Signature const&, Node);
  void InferBestSignature(Signature&, Signature&, Spfds&);

 private:
  Ntk& ntk_;
  AreaResynParams ps_;
  AreaResynStats& st_;
  WindowManagerParams wm_params_;
  SimulaManagerParams sm_params_;
  windowing::window_manager_stats wm_stats_;
  windowing::SimulaManagerStats sm_stats_;
  SharedOptoWindow window_;
  LutSynth synth_;
  SpfdMeth spfds_;
  Spfds spfd_;

  ArrivalsTracker arrivals_;
  std::vector<std::vector<kitty::dynamic_truth_table>> proj_fns_;
};

template <class Ntk, typename StaticParams>
void AreaResynImpl<Ntk, StaticParams>::ReportNtkState() {
  std::vector<uint32_t> mffc_sizes_solo;
  std::vector<uint32_t> mffc_sizes_dual;
  ntk_.foreach_gate([&](auto const& n) {
    window_.manager.run(n);
    uint32_t mffc_size = window_.manager.num_mffc();
    if (ntk_.num_outputs(n)==1)
    {
      while (mffc_sizes_solo.size() <= mffc_size) mffc_sizes_solo.push_back(0);
      mffc_sizes_solo[mffc_size] += 1;
    }
    else
    {
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
template <class Ntk, typename StaticParams>
void AreaResynImpl<Ntk, StaticParams>::Run() {
  const auto nmax = ntk_.size();
  int cnt = 0;

  ntk_.foreach_gate([&](auto const& n) {
    if (ntk_.is_pi(n)) return true;
    if (cnt++ >= nmax) return false;
    if (st_.GetNumOptos() >= ps_.stop_after_nth) return false;
    if (SkipNode(n)) return true;
    if (ntk_.num_outputs(n) != 1) return true;


    mockturtle::call_with_stopwatch(st_.time_window_constr, [&]() { window_.manager.run(n); });
    mockturtle::call_with_stopwatch(st_.time_window_simula, [&]() { window_.simulator.run(window_.manager);});

    OptimizeSolo(n);
    return true;
  });

  ntk_.foreach_gate([&](auto const& n) {
    if (ntk_.is_pi(n)) return true;
    if (cnt++ >= nmax) return false;
    if (st_.GetNumOptos() >= ps_.stop_after_nth) return false;
    if (SkipNode(n)) return true;
    if (ntk_.num_outputs(n) != 2) return true;

    mockturtle::call_with_stopwatch(st_.time_window_constr, [&]() { window_.manager.run(n); });
    mockturtle::call_with_stopwatch(st_.time_window_simula, [&]() { window_.simulator.run(window_.manager);});

    OptimizeDual(n);
    return true;
  });
}

template <class Ntk, typename StaticParams>
bool AreaResynImpl<Ntk, StaticParams>::SkipNode(Node n) const {
  if (ntk_.fanout_size(n) <= 0) return true;
  if (ntk_.is_pi(n) || ntk_.is_constant(n)) return true;
  if (ntk_.is_dead(n)) return true;
  return false;
}

template <class Ntk, typename StaticParams>
void AreaResynImpl<Ntk, StaticParams>::OptimizeSolo(Node n) {
  spfds_.Init();
  if (ps_.try_0resub_solo && mockturtle::call_with_stopwatch(st_.time_0resub_solo,
                                      [&]() { return Try0ResubSolo(n); })) {
    st_.Record0ResubSolo();
    return;
  }
  if (ps_.try_0rexub_solo && mockturtle::call_with_stopwatch(st_.time_0rexub_solo,
                                      [&]() { return Try0RexubSolo(n); })) {
    st_.Record0RexubSolo();
    return;
  }
  if (ps_.try_1resub_solo && mockturtle::call_with_stopwatch(st_.time_1resub_solo,
                                      [&]() { return Try1ResubSolo(n); })) {
    st_.Record1ResubSolo();
    return;
  }
  if (ps_.try_1rexub_solo && mockturtle::call_with_stopwatch(st_.time_1rexub_solo,
                                      [&]() { return Try1RexubSolo(n); })) {
    st_.Record1RexubSolo();
    return;
  }
  if (ps_.try_2resub_solo && mockturtle::call_with_stopwatch(st_.time_2resub_solo,
                                      [&]() { return Try2ResubSolo(n); })) {
    st_.Record2ResubSolo();
    return;
  }
  if (ps_.try_merger_solo && mockturtle::call_with_stopwatch(st_.time_merger_solo,
                                      [&]() { return TryMergerSolo(n); })) {
    st_.RecordMergerSolo();
    return;
  }
}

template <class Ntk, typename StaticParams>
bool AreaResynImpl<Ntk, StaticParams>::Try0ResubSolo(Node n) {
  std::optional<std::pair<Signal, bool>> resub;
  auto const& onset = window_.simulator.Get(ntk_.make_signal(n));
  auto const& offset = ~onset;
  window_.manager.foreach_divisor([&](auto const& f, auto) {
    auto const& divisor = window_.simulator.Get(f);
    if (kitty::equal(onset, divisor)) {
      if (!resub || resub->second) resub = std::pair{f, false};
    } else if (!resub && (ntk_.fanin_size(f.index) > 1) &&
               kitty::equal(offset, divisor))
      resub = std::pair{f, true};
  });
  if (resub) {
    if (resub->second)
      ntk_.substitute_node(n, ntk_.create_not(resub->first));
    else
      ntk_.substitute_node(n, resub->first);
    return true;
  }
  return false;
}

template <class Ntk, typename StaticParams>
bool AreaResynImpl<Ntk, StaticParams>::Try0RexubSolo(Node n) {
  // iterate one new node at a time
  size_t min_num_needed = 0;
  Signature target = window_.simulator.Get(ntk_.make_signal(n));
  std::vector<Signal> fanins;
  spfd_.init(target);

  uint32_t best_cost = std::numeric_limits<uint32_t>::max();
  std::optional<Node> best_div;
  window_.manager.foreach_divisor([&](auto const& fd, auto) {
    auto nd = ntk_.get_node(fd);
    if (!window_.manager.is_input(nd) && (ntk_.num_outputs(nd) == 1) &&
        (ntk_.fanin_size(nd) <= StaticParams::kMaxInpsDual)) {
      uint32_t const cost = spfd_.evaluate(window_.simulator.Get(fd));
      if (cost < best_cost) {
        best_cost = cost;
        best_div = nd;
      }
    }
  });
  if (best_div) {
    auto const fnew = AddOutput(*best_div, target);
    if (fnew) {
      window_.simulator.compute( window_.manager, (*fnew).index );
      bool cec1 = kitty::equal(window_.simulator.Get(Signal{(*fnew).index, 0}), window_.simulator.Get( ntk_.make_signal(*best_div) ));
      bool cec2 = kitty::equal(window_.simulator.Get(Signal{(*fnew).index, 1}), target );
      
      if (cec1 && cec2)
      {
        ntk_.substitute_node(*best_div, Signal{(*fnew).index, 0});
        ntk_.substitute_node(n, Signal{(*fnew).index, 1});

        return true;
      }
      else
      {
        std::cout << "Failed cec 0 Rexub " << std::endl;
        ntk_.take_out_node((*fnew).index);
      }
    }
  }
  return false;
}

template <class Ntk, typename StaticParams>
bool AreaResynImpl<Ntk, StaticParams>::Try1ResubSolo(Node n) {
  if (window_.manager.num_mffc() <= 1) return false;
  spfds_.Init();
  auto cut = spfds_.GetNext(StaticParams::kMaxInpsSolo);
  if (cut) {
    std::vector<Signal> children;
    for (auto f : *cut) children.push_back(f);

    auto tt = Interpolate(children, window_.simulator.Get(ntk_.make_signal(n)));
    auto fnew = ntk_.create_node(children, tt);
    window_.simulator.compute(window_.manager, fnew.index);
    if (kitty::equal(window_.simulator.Get(fnew), window_.simulator.Get(ntk_.make_signal(n))))
    {
      ntk_.substitute_node(n, fnew);
      return true;
    }
    else
    {
      ntk_.take_out_node(fnew.index);
      std::cout << "Failed cec 1 Resub " << std::endl;
    }
  }
  return false;
}

template <class Ntk, typename StaticParams>
bool AreaResynImpl<Ntk, StaticParams>::Try2ResubSolo(Node n) {
  if (window_.manager.num_mffc() <= 2) return false;
  spfd_.init(window_.simulator.Get(ntk_.make_signal(n)));
  std::vector<Signal> cutT;
  while ((cutT.size() < (StaticParams::kMaxInpsSolo-1))&&!spfd_.is_covered())
  {
    std::optional<Signal> d;
    uint32_t best_cost = std::numeric_limits<uint32_t>::max();
    window_.manager.foreach_divisor([&](auto const& fd, auto){
      auto nd = ntk_.get_node(fd);
      uint32_t const cost = spfd_.evaluate(window_.simulator.Get(fd));
      if (cost<best_cost)
      {
        best_cost = cost;
        d=fd;
      }
    });
    if(d)
    {
      cutT.push_back(*d);
      spfd_.update(window_.simulator.Get(*d));
    }
    else
      return false;
  }
  if (spfd_.is_covered()) return false;
  Signature funcB;
  Signature careB;
  InferBestSignature(funcB, careB, spfd_);

  spfd_.init(funcB, careB);
  std::vector<Signal> cutB;
  while (cutB.size() < StaticParams::kMaxInpsSolo &&!spfd_.is_covered())
  {
    std::optional<Signal> d;
    uint32_t best_cost = std::numeric_limits<uint32_t>::max();
    window_.manager.foreach_divisor([&](auto const& fd, auto){
      auto nd = ntk_.get_node(fd);
      uint32_t const cost = spfd_.evaluate(window_.simulator.Get(fd));
      if (cost<best_cost)
      {
        best_cost = cost;
        d=fd;
      }
    });
    if(d)
    {
      cutB.push_back(*d);
      spfd_.update(window_.simulator.Get(*d));
    }
    else
      return false;
  }
  if (spfd_.is_covered())
  {
    auto ttB = Interpolate(cutB, funcB, careB);
    auto fB = ntk_.create_node(cutB, ttB);
    window_.simulator.compute(window_.manager, fB.index);
    cutT.push_back(fB);
    auto ttT = Interpolate(cutT, window_.simulator.Get(ntk_.make_signal(n)));
    auto fT = ntk_.create_node(cutT, ttT);
    window_.simulator.compute(window_.manager, fT.index);
    if (kitty::equal(window_.simulator.Get(fT), window_.simulator.Get(ntk_.make_signal(n))))
    {
      ntk_.substitute_node(n, fT);
      return true;
    }
    else
    {
      std::cout << "Failed cec 2 Resub " << std::endl;
      ntk_.take_out_node(fT.index);
    }
  }

  return false;
}

template <class Ntk, typename StaticParams>
auto AreaResynImpl<Ntk, StaticParams>::AddOutput(
    Node n, Signature const& target) -> std::optional<Signal> {

  auto children = ntk_.get_children(n);

  spfd_.init(target);
  for (auto f : children) spfd_.update(window_.simulator.Get(f));

  if(children.size() < StaticParams::kMaxInpsDual)
  {
    ntk_.incr_trav_id();
    window_.manager.MarkTfoInside(n, ntk_.trav_id());
    auto fn = ntk_.make_signal(n);
    while ((children.size() < StaticParams::kMaxInpsDual)&&!spfd_.is_covered()) {
      std::optional<Signal> best_div;
      uint32_t best_cost = std::numeric_limits<uint32_t>::max();
      window_.manager.foreach_divisor([&](auto const& fd, auto){
        auto nd = ntk_.get_node(fd);
        if ((ntk_.visited(nd)!=ntk_.trav_id())&&((arrivals_.get_time(fd)<arrivals_.get_time(fn)) ))//|| !ExistsPath(n,nd)))
        {
          uint32_t const cost = spfd_.evaluate(window_.simulator.Get(fd));
          if (cost<best_cost)
          {
            best_cost = cost;
            best_div = fd;
          }
        }
      });
      if ( best_div )
      {
        children.push_back(*best_div);
        spfd_.update(window_.simulator.Get(*best_div) );
      }
      else 
        return std::nullopt;
    }
  }

  if (spfd_.is_covered()) {
    auto tt1 = kitty::extend_to(ntk_.node_function(n, 0), children.size());
    auto tt2 = Interpolate(children, target);
    auto fnew = ntk_.create_node(children, {tt1, tt2});

    return fnew;
  }
  return std::nullopt;
}

template <class Ntk, typename StaticParams>
auto AreaResynImpl<Ntk, StaticParams>::AddOutput(
    Node n, Signature const& func, Signature const& care) -> std::optional<Signal> {

  auto children = ntk_.get_children(n);

  spfd_.init(func, care);
  for (auto f : children) spfd_.update(window_.simulator.Get(f));

  if(children.size() < StaticParams::kMaxInpsDual)
  {
    ntk_.incr_trav_id();
    window_.manager.MarkTfoInside(n, ntk_.trav_id());
    auto fn = ntk_.make_signal(n);
    while ((children.size() < StaticParams::kMaxInpsDual)&&!spfd_.is_covered()) {
      std::optional<Signal> best_div;
      uint32_t best_cost = std::numeric_limits<uint32_t>::max();
      window_.manager.foreach_divisor([&](auto const& fd, auto){
        auto nd = ntk_.get_node(fd);
        if ((ntk_.visited(nd)!=ntk_.trav_id())&&((arrivals_.get_time(fd)<arrivals_.get_time(fn)) ))
        {
          uint32_t const cost = spfd_.evaluate(window_.simulator.Get(fd));
          if (cost<best_cost)
          {
            best_cost = cost;
            best_div = fd;
          }
        }
      });
      if ( best_div )
      {
        children.push_back(*best_div);
        spfd_.update(window_.simulator.Get(*best_div) );
      }
      else 
        return std::nullopt;
    }
  }

  if (spfd_.is_covered()) {
    auto tt1 = kitty::extend_to(ntk_.node_function(n, 0), children.size());
    auto tt2 = Interpolate(children, func, care);
    auto fnew = ntk_.create_node(children, {tt1, tt2});

    return fnew;
  }
  return std::nullopt;
}

template <class Ntk, typename StaticParams>
bool AreaResynImpl<Ntk, StaticParams>::TryMergerSolo(Node n) {
  // Collect inputs of n
  std::set<Signal> base_inputs;
  ntk_.foreach_fanin(n, [&](auto const& fi, auto) { base_inputs.insert(fi); });

  std::optional<Node> best_node;
  int best_cost = std::numeric_limits<int>::max();

  // Explore fanouts of fanins
  ntk_.foreach_fanin(n, [&](auto const& fi, auto) {
    ntk_.foreach_fanout(fi.index, [&](auto const& no) {
      if (n == no) return;
      if (ntk_.num_outputs(no) != 1u) return;

      // Compute combined inputs size
      std::set<Signal> merged_inputs = base_inputs;
      ntk_.foreach_fanin(
          no, [&](auto const& foi, auto) { merged_inputs.insert(foi); });

      if (!ExistsPath(n, no) && ((int)merged_inputs.size() < best_cost)) {
        {
          best_node = no;
          best_cost = (int)merged_inputs.size();
        }
      }
    });
  });

  if (best_node) {
    Merge(n, *best_node);
    return true;
  }
  return false;
}

template <class Ntk, typename StaticParams>
bool AreaResynImpl<Ntk, StaticParams>::ExistsPath(Node n1, Node n2) const {
  if (n1 == n2) return true;
  Node n_upper = n1;
  Node n_lower = n2;
  Signal f_upper = ntk_.make_signal(n_upper);
  Signal f_lower = ntk_.make_signal(n_lower);
  double t_upper = arrivals_.get_time(f_upper);
  double t_lower = arrivals_.get_time(f_lower);

  if (t_upper < t_lower) {
    std::swap(n_upper, n_lower);
    std::swap(f_upper, f_lower);
  }

  ntk_.incr_trav_id();
  bool found = false;
  std::function<void(Node)> DepthFirstSearch = [&](Node u) {
    if (found) return;
    if (u == n_lower) {
      found = true;
      return;
    }
    ntk_.foreach_fanin(u, [&](auto fi, auto) {
      auto v = ntk_.get_node(fi);
      if (ntk_.visited(v) != ntk_.trav_id() &&
          (arrivals_.get_time(fi) >= (arrivals_.get_time(f_lower) * 0.9))) {
        ntk_.set_visited(v, ntk_.trav_id());
        DepthFirstSearch(v);
      }
    });
  };
  DepthFirstSearch(n_upper);
  return found;
}

template <class Ntk, typename StaticParams>
void AreaResynImpl<Ntk, StaticParams>::Merge(Node n1, Node n2) {
  if (ntk_.is_pi(n1) || ntk_.is_pi(n2)) return;

  std::unordered_map<uint64_t, std::pair<int, int>> map;
  std::set<Signal> inps;
  ntk_.foreach_fanin(n1, [&](auto const& fi, auto ii) { inps.insert(fi);});
  ntk_.foreach_fanin(n2, [&](auto const& fi, auto ii) { inps.insert(fi);});

  if (inps.size() > StaticParams::kMaxInpsDual) return;
  int i = 0;
  std::vector<Signal> children;
  for (auto const& f : inps) {
    children.push_back(f);
    map[f.data] = std::pair{inps.size(), i++};
  }
  std::vector<kitty::dynamic_truth_table const*> sim_ptrs1;
  std::vector<kitty::dynamic_truth_table const*> sim_ptrs2;

  ntk_.foreach_fanin(n1, [&](auto const& fi, auto ii) {
    auto const p = map[fi.data];
    sim_ptrs1.push_back(&proj_fns_[p.first][p.second]);
  });
  auto funcs1 = ntk_.compute(n1, sim_ptrs1);

  ntk_.foreach_fanin(n2, [&](auto const& fi, auto ii) {
    auto const p = map[fi.data];
    sim_ptrs2.push_back(&proj_fns_[p.first][p.second]);
  });
  auto funcs2 = ntk_.compute(n2, sim_ptrs2);

  std::vector<kitty::dynamic_truth_table> tts = {funcs1[0], funcs2[0]};

  auto fnew = ntk_.create_node(children, tts);

  ntk_.substitute_node(n1, Signal{fnew.index, 0});
  ntk_.substitute_node(n2, Signal{fnew.index, 1});
}

template <class Ntk, typename StaticParams>
kitty::dynamic_truth_table AreaResynImpl<Ntk, StaticParams>::Interpolate(
    std::vector<Signal> const& children, Signature const& target) {
  uint32_t const num_vars = children.size();
  kitty::dynamic_truth_table tt(num_vars);
  uint32_t const num_minterms = 1u << num_vars;
  Signature mask;
  for (uint32_t m = 0; m < num_minterms; ++m) {
    mask |= ~mask;
    for (uint32_t v = 0; v < num_vars; ++v) {
      if (((m >> v) & 0x1) > 0)
        mask &= window_.simulator.Get(children[v]);
      else
        mask &= ~window_.simulator.Get(children[v]);
    }
    uint32_t c1 = kitty::count_ones(mask & target);
    uint32_t c0 = kitty::count_ones(mask&~target);
    if (c0==0&&c1>0)
      kitty::set_bit(tt, m);
    else if (c1==0)
      kitty::clear_bit(tt, m);
    else
    {
      std::cout << "Warning: don't know how to interpolate" << std::endl;
    }
  }
  return tt;
}

template <class Ntk, typename StaticParams>
kitty::dynamic_truth_table AreaResynImpl<Ntk, StaticParams>::Interpolate(
    std::vector<Signal> const& children, Signature const& target, Signature const& care) {

  uint32_t const num_vars = children.size();
  kitty::dynamic_truth_table tt(num_vars);
  uint32_t const num_minterms = 1u << num_vars;
  Signature mask;
  for (uint32_t m = 0; m < num_minterms; ++m) {
    mask |= ~mask;
    for (uint32_t v = 0; v < num_vars; ++v) {
      if (((m >> v) & 0x1) > 0)
        mask &= window_.simulator.Get(children[v]);
      else
        mask &= ~window_.simulator.Get(children[v]);
    }
    if (kitty::count_ones(mask & target & care) > kitty::count_ones(mask & care & ~target))
      kitty::set_bit(tt, m);
    else
      kitty::clear_bit(tt, m);
  }
  return tt;
}

template <class Ntk, typename StaticParams>
bool AreaResynImpl<Ntk, StaticParams>::Try1RexubSolo(Node n) {
  return false;
  if (window_.manager.num_mffc() <= 1) return false;
  Signature func = window_.simulator.Get(ntk_.make_signal(n));
  Signature care = ~func.construct();

  std::vector<Node> added; 
  std::vector<std::pair<Node, Signal>> todos;
  std::vector<Signal> new_fanin;

  ntk_.incr_trav_id();
  spfd_.init(func);
  while ( (new_fanin.size()<StaticParams::kMaxInpsSolo)&&!(spfd_.is_covered()) )
  {
    std::optional<Signal> best_asis;
    std::optional<Node> best_dupl;
    uint32_t best_cost_asis = std::numeric_limits<uint32_t>::max();
    uint32_t best_cost_dupl = std::numeric_limits<uint32_t>::max();
    window_.manager.foreach_divisor([&](auto const& fd, auto){
      auto const nd = ntk_.get_node(fd);
      if (ntk_.visited(nd)!=ntk_.trav_id())
      {
        uint32_t const cost = spfd_.evaluate(window_.simulator.Get(fd));
        if (cost<best_cost_asis)
        {
          best_cost_asis=cost;
          best_asis = fd;
        }
        if ((!window_.manager.is_input(nd))&&(ntk_.num_outputs(nd)==1)&&(cost<best_cost_dupl)&&(ntk_.fanin_size(nd)<=StaticParams::kMaxInpsDual))
        {
          best_cost_dupl=cost;
          best_dupl = nd;
        }
      }
    });
    if (best_dupl)
    {
      auto fnew = AddOutput(*best_dupl, func, care);
      if (fnew)
      {
        Node nnew = ntk_.get_node(*fnew);
        window_.simulator.compute(window_.manager, nnew);
        auto cost = spfd_.evaluate(window_.simulator.Get(Signal{nnew, 1}));
        if (cost<best_cost_asis)
        {
          added.push_back(nnew);
          todos.push_back(std::pair{*best_dupl, Signal{nnew, 0}});
          new_fanin.push_back(Signal{nnew, 1});
          spfd_.update(window_.simulator.Get(Signal{nnew, 1}));
          ntk_.set_visited(*best_dupl, ntk_.trav_id());
        }
        else if (best_asis)
        {
          if (ntk_.fanout_size(nnew)==0) ntk_.take_out_node(nnew);
  
          new_fanin.push_back(*best_asis);
          spfd_.update(window_.simulator.Get(*best_asis));
        }
      }
    }
    else if (best_asis)
    {
      new_fanin.push_back(*best_asis);
      spfd_.update(window_.simulator.Get(*best_asis));
    }
    InferBestSignature( func, care, spfd_ );
  }
  
  if(spfd_.is_covered())
  {
      auto tt = Interpolate(new_fanin, window_.simulator.Get(ntk_.make_signal(n)) );
      auto fnew = ntk_.create_node( new_fanin, tt );
      window_.simulator.compute(window_.manager, fnew.index);
      if(!kitty::equal(window_.simulator.Get(ntk_.make_signal(n)), window_.simulator.Get(fnew))){
        std::cout << "Failed cec 1 Rexub " << std::endl;
        ntk_.take_out_node(fnew.index);
        return false;
      }
      ntk_.substitute_node(n, fnew);
      for (auto const p : todos)
        ntk_.substitute_node( p.first, p.second );
      return true;
  }
  for (auto nd : added){ 
    if (ntk_.fanout_size(nd)==0) ntk_.take_out_node(nd);
  }
  return false;
}

template <class Ntk, typename StaticParams>
kitty::dynamic_truth_table AreaResynImpl<Ntk, StaticParams>::Approximate(Signature const& target, Signature const& care, Node n)
{
  Signature temp;
  uint32_t num_vars = ntk_.fanin_size(n);
  kitty::dynamic_truth_table tt(num_vars);
  uint32_t num_minterms = 1u << num_vars;
  auto const children = ntk_.get_children(n);
  for (uint32_t m=0u; m<num_minterms; ++m)
  {
    temp |= ~temp;
    for (uint32_t v=0; v<num_vars; ++v)
      temp &= ((m>>v)&0x1) ? window_.simulator.Get(children[v]) : ~window_.simulator.Get(children[v]);
    
    uint32_t const num_1s = kitty::count_ones(temp&target&care);
    uint32_t const num_0s = kitty::count_ones(temp&(~target)&care);
    if ( num_1s > num_0s )
      kitty::set_bit(tt,m);
    else
      kitty::clear_bit(tt,m);
  }
  return tt;
}

template <class Ntk, typename StaticParams>
void AreaResynImpl<Ntk, StaticParams>::InferBestSignature(Signature&func, Signature&care, Spfds& spfd)
{
  func ^= func;
  care ^= care;
  for ( int m=0; m<spfd.num_masks_; ++m )
  {
    uint32_t const num1=kitty::count_ones(spfd.masks_[m]&spfd.func_[1]);
    uint32_t const num0=kitty::count_ones(spfd.masks_[m]);
    if ((num1*num0)!=0) care |= spfd.masks_[m];

    if (kitty::count_ones(spfd.masks_[m]&spfd.func_[1])>kitty::count_ones(spfd.masks_[m]&spfd.func_[0]))
      func |= spfd.masks_[m] & spfd.func_[1];
    else
      func |= spfd.masks_[m] & spfd.func_[0];
  }
}

template <class Ntk, typename StaticParams>
void AreaResynImpl<Ntk, StaticParams>::OptimizeDual(Node n) {
  if (ps_.try_1resub_dual && mockturtle::call_with_stopwatch(st_.time_1resub_dual,
                                      [&]() { return Try1ResubDual(n); })) {
    st_.Record1ResubDual();
    return;
  }
}

template <class Ntk, typename StaticParams>
bool AreaResynImpl<Ntk, StaticParams>::Try1ResubDual(Node n) {
  if (window_.manager.num_mffc() <= 1) return false;
  spfds_.Init();
  auto cut = spfds_.GetNext(StaticParams::kMaxInpsDual);
  if (cut) {
    std::vector<Signal> children;
    for (auto f : *cut) {
      children.push_back(f);
    }

    auto tt0 = Interpolate(children, window_.simulator.Get(Signal{n,0}));
    auto tt1 = Interpolate(children, window_.simulator.Get(Signal{n,1}));
    auto fnew = ntk_.create_node(children, {tt0, tt1});
    auto nnew = ntk_.get_node(fnew);
    window_.simulator.compute(window_.manager, nnew);
    bool cec1 = kitty::equal(window_.simulator.Get(Signal{n,0}), window_.simulator.Get(Signal{nnew,0}));
    bool cec2 = kitty::equal(window_.simulator.Get(Signal{n,1}), window_.simulator.Get(Signal{nnew,1}));

    if (cec1 && cec2)
    {
      std::vector<Signal> new_signals = { Signal{nnew, 0}, Signal{nnew, 1} };
      ntk_.substitute_node(n, new_signals);
      return true;
    }
    {
      std::cout << "Failed cec 1 Resub dual " << std::endl;
      ntk_.take_out_node(nnew);
    }
  }
  return false;
}


}  // namespace detail

// Public API
template <class Ntk, typename StaticParams = DefaultAreaResynStaticParams>
void AreaResyn(Ntk& ntk, AreaResynParams ps = {},
               AreaResynStats* pst = nullptr) {
  AreaResynStats st;
  detail::AreaResynImpl<Ntk, StaticParams> p(ntk, ps, st);
  p.ReportNtkState();
  p.Run();
  p.ReportNtkState();
  if (pst) *pst = st;
}

}  // namespace array
}  // namespace opto
}  // namespace rinox

#endif  // RINOX_ALGORITHMS_ARRAY_AREA_ORCHESTRA_HPP_
